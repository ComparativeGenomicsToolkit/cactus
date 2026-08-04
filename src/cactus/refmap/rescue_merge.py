#!/usr/bin/env python3
"""
Gap-fill rescue: fill query regions that minigraph left unaligned (after its normal
overlap filter) with direct assembly->reference alignments (minimap2 / cactus-refmap),
lifted into graph node coordinates so they merge into the graphmap PAF cleanly.

Pipeline per gap:
  - `paffy to_bed -f` finds the query coverage gaps in the graphmap PAF.
  - the uniqueness gate (primary cover, no query-overlapping primaries, no strong secondary)
    picks the reference records that map a gap cleanly -- run on the ORIGINAL refmap, which
    still carries the tp:A: tags.
  - `paffy shatter` breaks the selected records into gapless pieces (query-len == target-len).
  - each gapless piece is lifted to node coordinates against the reference->node table (built
    from the reference's own graphmap records), splitting at backbone-node boundaries by simple
    1:1 arithmetic -- an inverted piece becomes a reverse traversal of the backbone node.
The lifted, node-targeted fills are appended to the graphmap PAF (its records pass through
unchanged), and split/align consume the result normally.

Query names must match between the minigraph PAF and the reference PAF (id=<event>|<contig>).
"""
import argparse, sys, gzip, subprocess
from collections import defaultdict


def _open(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def read_paf(path):
    """query -> list of (qs, qe, tp, target, raw_line).  tp in {'P','S','I'}."""
    by_q = defaultdict(list)
    with _open(path) as f:
        for line in f:
            if not line.strip() or line[0] == '@':
                continue
            t = line.rstrip("\n").split("\t")
            if len(t) < 12:
                continue
            q, qs, qe, target = t[0], int(t[2]), int(t[3]), t[5]
            tp = 'P'
            for tag in t[12:]:
                if tag.startswith("tp:A:"):
                    tp = tag[5]
                    break
            by_q[q].append((qs, qe, tp, target, line))
    return by_q


def build_ref_node_table(graphmap_paf, ref_prefix):
    """From a graphmap PAF, read the reference's own records (query starts with ref_prefix, e.g.
    'id=ref|') to get the reference->node mapping straight from the pipeline.  Returns
    {ref_target: sorted [(ref_qs, ref_qe, node, node_len, node_ts, node_te)]}.  The reference is
    the rank-0 backbone, so it maps forward and gaplessly to its nodes (ref pos p -> node
    node_ts + (p - ref_qs))."""
    table = defaultdict(list)
    with _open(graphmap_paf) as f:
        for line in f:
            t = line.rstrip("\n").split("\t")
            if len(t) < 12 or not t[0].startswith(ref_prefix) or t[4] != '+':
                continue
            table[t[0]].append((int(t[2]), int(t[3]), t[5], int(t[6]), int(t[7]), int(t[8])))
    for k in table:
        table[k].sort()
    return table


def cluster_by_overlap(intervals, min_overlap):
    clusters, cur, cur_end = [], [], None
    for qs, qe, payload in sorted(intervals, key=lambda x: x[0]):
        if cur and qs < cur_end and (min(qe, cur_end) - qs) >= min_overlap:
            cur.append(payload)
            cur_end = max(cur_end, qe)
        else:
            if cur:
                clusters.append(cur)
            cur, cur_end = [payload], qe
    if cur:
        clusters.append(cur)
    return clusters


def covered_fraction(intervals, cs, ce):
    span = ce - cs
    if span <= 0:
        return 0.0
    clipped = sorted((max(qs, cs), min(qe, ce)) for qs, qe in intervals if min(qe, ce) > max(qs, cs))
    covered, s, e = 0, None, None
    for a, b in clipped:
        if e is None:
            s, e = a, b
        elif a <= e:
            e = max(e, b)
        else:
            covered += e - s
            s, e = a, b
    if e is not None:
        covered += e - s
    return covered / span


def rescue_records(ref_recs, cs, ce, cover_frac, secondary_frac, min_overlap):
    """Return the (assembly->reference) refmap lines to fill span [cs,ce], or None if the reference
    doesn't map it cleanly (needs a primary cover >= cover_frac, no primaries overlapping each
    other in query, and no strong secondary)."""
    prim, sec = [], []
    for qs, qe, tp, target, line in ref_recs:
        if min(qe, ce) <= max(qs, cs):
            continue
        (prim if tp in ('P', 'I') else sec).append((qs, qe, line))
    if not prim:
        return None
    if any(len(c) > 1 for c in cluster_by_overlap([(qs, qe, line) for qs, qe, line in prim], min_overlap)):
        return None                                   # reference maps the span to >1 place
    span = ce - cs
    if any((min(qe, ce) - max(qs, cs)) >= secondary_frac * span for qs, qe, _ in sec):
        return None                                   # a strong secondary => ambiguous
    if covered_fraction([(qs, qe) for qs, qe, _ in prim], cs, ce) < cover_frac:
        return None                                   # reference doesn't cover the span
    return [line for _, _, line in prim]


def to_bed_gaps(paf_path, min_gap, paffy="paffy", assembly_fa=None):
    """{query: [(gap_start, gap_end), ...]} -- query intervals with zero alignment coverage
    (length >= min_gap), via `paffy to_bed -f`."""
    cmd = [paffy, "to_bed", "-i", paf_path, "-f", "-m", str(min_gap)]
    if assembly_fa:
        cmd += ["-q", assembly_fa]
    out = subprocess.run(cmd, capture_output=True, text=True, check=True).stdout
    gaps = defaultdict(list)
    for line in out.splitlines():
        t = line.split()
        if len(t) >= 3:
            gaps[t[0]].append((int(t[1]), int(t[2])))
    return gaps


def lift_gapless_piece(t, table):
    """Lift one GAPLESS assembly->ref PAF record (query-len == target-len, no indels) to
    assembly->node record(s), splitting at backbone-node boundaries.  Linear 1:1 arithmetic;
    an inverted piece stays reverse (a reverse traversal of the node).

    The reference->node table records are NOT guaranteed gapless: minigraph's graphmap can align the
    reference to a node across an indel, so a record's ref-span (rqe-rqs) can exceed its node-span
    (nte-nts).  A naive 1:1 lift then runs the node coordinate off the end of the node
    (n_s = nts + (ov_s-rqs) > nlen), which cactus_consolidated rejects ("Paf target start
    coordinates are invalid").  So clamp the lifted node interval to the node span and clip the query
    to keep query-len == node-len -- fills past the node boundary are simply dropped."""
    if len(t) < 9 or t[5] not in table:
        return []
    q, qlen, qs, qe, strand = t[0], t[1], int(t[2]), int(t[3]), t[4]
    ts, te = int(t[7]), int(t[8])
    mapq = t[11] if len(t) > 11 else "60"
    out = []
    for (rqs, rqe, node, nlen, nts, nte) in table[t[5]]:
        ov_s, ov_e = max(ts, rqs), min(te, rqe)
        if ov_e <= ov_s:
            continue
        n_s = nts + (ov_s - rqs)
        if n_s >= nte:                       # overlap starts past this node (ref-span > node-span)
            continue
        seg = min(ov_e - ov_s, nte - n_s)    # clamp to the node span so we never overrun the node
        n_e = n_s + seg
        if strand == '+':
            q_s = qs + (ov_s - ts)
            q_e = q_s + seg
        else:
            q_e = qe - (ov_s - ts)
            q_s = q_e - seg
        out.append("\t".join([q, qlen, str(q_s), str(q_e), strand, node, str(nlen),
                              str(n_s), str(n_e), str(seg), str(seg), mapq,
                              "tp:A:P", "cg:Z:{}M".format(seg)]) + "\n")
    return out


def clip_piece_to_gaps(t, gaps):
    """Clip a gapless node-targeted PAF record (from lift_gapless_piece) to the query's gap
    intervals, so fills never overhang into minigraph-covered flanks (an overhang would give the
    contig both a forward and a reverse mapping to the same node -- a palindrome rgfa-split axes).
    Linear 1:1; handles reverse (query up => node down)."""
    qs, qe, strand = int(t[2]), int(t[3]), t[4]
    ns, ne = int(t[7]), int(t[8])
    out = []
    for gs, ge in gaps:
        iq_s, iq_e = max(qs, gs), min(qe, ge)
        if iq_e <= iq_s:
            continue
        if strand == '+':
            n_s, n_e = ns + (iq_s - qs), ns + (iq_e - qs)
        else:
            n_s, n_e = ne - (iq_e - qs), ne - (iq_s - qs)
        seg = iq_e - iq_s
        out.append("\t".join([t[0], t[1], str(iq_s), str(iq_e), strand, t[5], t[6],
                              str(n_s), str(n_e), str(seg), str(seg), t[11],
                              "tp:A:P", "cg:Z:{}M".format(seg)]) + "\n")
    return out


def gap_fill(minigraph_paf, rm_by_q, ref_table, out_paf, min_gap=1000, cover_frac=0.5,
             secondary_frac=0.5, paffy="paffy", assembly_fa=None):
    """Write out_paf = the minigraph PAF unchanged + node-targeted reference records that fill its
    query coverage gaps.  rm_by_q is the ORIGINAL assembly->reference refmap (read_paf); ref_table
    is the reference->node table (build_ref_node_table).  Returns a stats dict."""
    gaps = to_bed_gaps(minigraph_paf, min_gap, paffy, assembly_fa)
    st = dict(gaps=0, filled=0, fill_records=0, unfilled=0)
    # gate on the original refmap (tp:A: tags intact) -> the reference lines that cleanly fill gaps
    selected = []
    for q, ivs in gaps.items():
        for gs, ge in ivs:
            st["gaps"] += 1
            resc = rescue_records(rm_by_q.get(q, []), gs, ge, cover_frac, secondary_frac, min_gap)
            if resc:
                selected.extend(resc)
                st["filled"] += 1
            else:
                st["unfilled"] += 1
    # shatter the selected records into gapless pieces, then lift each to node coordinates
    fills = []
    if selected:
        sel_path = out_paf + ".selected.paf"
        with open(sel_path, "w") as f:
            f.writelines(selected)
        shattered = subprocess.run([paffy, "shatter", "-i", sel_path],
                                   capture_output=True, text=True, check=True).stdout
        for line in shattered.splitlines():
            for piece in lift_gapless_piece(line.split("\t"), ref_table):
                pt = piece.rstrip("\n").split("\t")
                fills.extend(clip_piece_to_gaps(pt, gaps.get(pt[0], [])))   # keep fills inside the gaps only
    # final safety net: never emit a fill whose node interval is out of bounds or whose query/node
    # spans disagree -- cactus_consolidated rejects such records and the whole run dies hours later.
    # Should be a no-op given lift_gapless_piece clamps to the node, but 200k fills * one bad record
    # is not worth the risk.
    valid, dropped = [], 0
    for line in fills:
        c = line.rstrip("\n").split("\t")
        nlen, n_s, n_e = int(c[6]), int(c[7]), int(c[8])
        if 0 <= n_s < n_e <= nlen and (int(c[3]) - int(c[2])) == (n_e - n_s):
            valid.append(line)
        else:
            dropped += 1
    fills = valid
    st["fill_records"] = len(fills)
    st["fill_dropped"] = dropped
    with open(out_paf, "w") as out:
        with _open(minigraph_paf) as f:               # minigraph records pass through unchanged
            for line in f:
                out.write(line if line.endswith("\n") else line + "\n")
        for line in fills:
            out.write(line)
    return st


def gap_fill_pafs(minigraph_paf, refmap_paf, out_paf, ref_prefix, min_gap=1000, cover_frac=0.5,
                  secondary_frac=0.5, paffy="paffy", assembly_fa=None):
    """Convenience: build the ref->node table from the minigraph PAF and gap-fill one file."""
    table = build_ref_node_table(minigraph_paf, ref_prefix)
    return gap_fill(minigraph_paf, read_paf(refmap_paf), table, out_paf, min_gap, cover_frac,
                    secondary_frac, paffy, assembly_fa)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--graphmap", required=True, help="minigraph graphmap PAF (may be .gz)")
    ap.add_argument("--refmap", required=True, help="assembly->reference PAF from cactus-refmap (may be .gz)")
    ap.add_argument("-o", "--output", required=True, help="output merged PAF")
    ap.add_argument("--ref-prefix", required=True, help="reference query prefix, e.g. 'id=CHM13|'")
    ap.add_argument("--min-gap", type=int, default=1000)
    ap.add_argument("--cover-frac", type=float, default=0.5)
    ap.add_argument("--secondary-frac", type=float, default=0.5)
    ap.add_argument("--paffy", default="paffy", help="path to the paffy binary [paffy]")
    ap.add_argument("--assembly", default=None, help="query fasta, to also surface fully-missing contigs (paffy -q)")
    args = ap.parse_args()
    st = gap_fill_pafs(args.graphmap, args.refmap, args.output, args.ref_prefix, args.min_gap,
                       args.cover_frac, args.secondary_frac, args.paffy, args.assembly)
    sys.stderr.write("[rescue] {gaps} gaps: {filled} filled -> {fill_records} node records, "
                     "{unfilled} left unfilled\n".format(**st))


if __name__ == "__main__":
    main()
