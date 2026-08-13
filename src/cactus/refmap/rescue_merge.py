#!/usr/bin/env python3
"""
Gap-fill rescue: fill query regions that minigraph left unaligned (after its normal
overlap filter) with direct assembly->reference alignments (minimap2 / cactus-refmap),
lifted into graph node coordinates so they merge into the graphmap PAF cleanly.

Pipeline per gap:
  - find_query_gaps finds the query coverage gaps in the graphmap PAF (interval scan, not paffy to_bed).
  - the uniqueness gate (primary cover, no query-overlapping primaries, no strong secondary)
    picks the reference records that map a gap cleanly -- run on the ORIGINAL refmap, which
    still carries the tp:A: tags.
  - `paffy shatter` breaks the selected records into gapless pieces (query-len == target-len).
  - each gapless piece is (optionally) gated on the both-haplotype single-coverage confidence set
    -- dipcall's SD filter, rebuilt from the refmap (build_confident_intervals) -- then lifted to
    node coordinates against the reference->node table (built from the reference's own graphmap
    records), splitting at backbone-node boundaries by simple 1:1 arithmetic -- an inverted piece
    becomes a reverse traversal of the backbone node.
The lifted, node-targeted fills are appended to the graphmap PAF (its records pass through
unchanged), and split/align consume the result normally.

Query names must match between the minigraph PAF and the reference PAF (id=<event>|<contig>).
"""
import argparse, sys, gzip, subprocess, re, bisect
from collections import defaultdict

_CIGAR_RE = re.compile(r'(\d+)([MIDX=])')


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


def rescue_records(ref_recs, cs, ce, cover_frac, secondary_frac, min_overlap, min_mapq=0, min_block=0):
    """Return the (assembly->reference) refmap lines to fill span [cs,ce], or None if the reference
    doesn't map it cleanly (needs a primary cover >= cover_frac, no primaries overlapping each other
    WITHIN the gap, and no strong secondary).

    dipcall's samflt gate: a gap is only filled from source alignments that are primary, mapQ >=
    min_mapq, and >= min_block bp of alignment block (PAF col 11).  This is the mapQ/length filter
    the minigraph path already applies (filter_paf minMAPQ/minGAFBlockLength); on the minimap2 refmap
    it drops the short/low-mapQ paralog fragments that are ~99% of the SD collateral -- validated:
    refmap alignments <50kb are 99% in dipcall-ambiguous regions, >=50kb are ~2/3 real."""
    prim, sec = [], []
    for qs, qe, tp, target, line in ref_recs:
        if min(qe, ce) <= max(qs, cs):
            continue
        if min_mapq or min_block:
            t = line.rstrip("\n").split("\t")
            if len(t) < 12 or int(t[11]) < min_mapq or int(t[10]) < min_block:
                continue                       # samflt: drop low-mapQ / short source alignments
        (prim if tp in ('P', 'I') else sec).append((qs, qe, line))
    if not prim:
        return None
    # test "maps to >1 place" on each primary's portion WITHIN the gap, not its full extent: an
    # SD-mediated inversion aligns the query's inverted core reverse to the reference while its forward
    # flanks (the inverted-repeat arms) overlap that core OUTSIDE the gap -- clipping to [cs,ce] stops
    # that spurious overlap from rejecting a clean single-primary fill.  A genuine paralog duplicate
    # maps two primaries onto the SAME gap span, so it still overlaps in-gap and is still rejected.
    if any(len(c) > 1 for c in cluster_by_overlap(
            [(max(qs, cs), min(qe, ce), line) for qs, qe, line in prim], min_overlap)):
        return None                                   # reference maps the gap to >1 place
    span = ce - cs
    if any((min(qe, ce) - max(qs, cs)) >= secondary_frac * span for qs, qe, _ in sec):
        return None                                   # a strong secondary => ambiguous
    if covered_fraction([(qs, qe) for qs, qe, _ in prim], cs, ce) < cover_frac:
        return None                                   # reference doesn't cover the span
    return [line for _, _, line in prim]


def _covered_query_intervals(qs, qe, strand, cigar):
    """Query intervals actually ALIGNED by the cigar (maximal runs of M/=/X).  An insertion (I)
    consumes query without aligning it, so it breaks a run and reads as a coverage gap -- exactly the
    unaligned query sequence the rescue targets (e.g. an inverted middle minigraph couldn't place).
    D/N consume target only, so they leave query coverage contiguous.  Reverse strand consumes the
    query from qe downward, matching paffy's coverage."""
    ivs = []
    if strand != '-':
        q = run_start = qs
        for n, op in _CIGAR_RE.findall(cigar):
            n = int(n)
            if op == 'I':
                if q > run_start:
                    ivs.append((run_start, q))
                q += n
                run_start = q
            elif op != 'D':                        # M/=/X advance & cover; D leaves the query alone
                q += n
        if q > run_start:
            ivs.append((run_start, q))
    else:
        q = run_end = qe
        for n, op in _CIGAR_RE.findall(cigar):
            n = int(n)
            if op == 'I':
                if run_end > q:
                    ivs.append((q, run_end))
                q -= n
                run_end = q
            elif op != 'D':
                q -= n
        if run_end > q:
            ivs.append((q, run_end))
    return ivs


def find_query_gaps(paf_path, min_gap):
    """{query: [(gap_start, gap_end), ...]} -- maximal query intervals with ZERO alignment coverage
    and length >= min_gap.  Reproduces `paffy to_bed -f -m <min_gap>` (minus the -q fasta option the
    rescue never uses) but from the alignment CIGARs instead of a per-base coverage array: to_bed
    allocates 2 bytes for every base of every query and holds them all at once (~68 GiB on a 7-way,
    ~TBs at HPRC scale); this is O(records + insertions), not O(query bases).  Coverage is the union
    of aligned (M/=/X) query intervals, so an unaligned insertion reads as a gap -- the rescue's
    signal.  It is NOT the raw [qs,qe] span: a record can span a region it doesn't actually align."""
    covered = defaultdict(list)
    qlen = {}
    with _open(paf_path) as f:
        for line in f:
            if not line.strip() or line[0] == '@':
                continue
            t = line.rstrip("\n").split("\t")
            if len(t) < 12:
                continue
            q, ql, qs, qe, strand = t[0], int(t[1]), int(t[2]), int(t[3]), t[4]
            qlen[q] = ql
            cg = next((tag[5:] for tag in t[12:] if tag.startswith("cg:Z:")), None)
            if cg:
                covered[q].extend(_covered_query_intervals(qs, qe, strand, cg))
            else:
                covered[q].append((qs, qe))        # no cigar: treat the whole span as aligned
    gaps = defaultdict(list)
    for q, ivs in covered.items():
        ivs.sort()
        pos = 0                                    # end of the merged coverage so far
        for s, e in ivs:
            if s > pos:                            # [pos, s) is uncovered
                if s - pos >= min_gap:
                    gaps[q].append((pos, s))
                pos = e
            elif e > pos:
                pos = e
        if qlen[q] - pos >= min_gap:               # trailing gap after the last alignment
            gaps[q].append((pos, qlen[q]))
    return gaps


def find_flanked_gaps(paf_path, min_gap):
    """find_query_gaps + each gap's immediate flank lengths: {query: [(gs, ge, lflank, rflank)]}.
    lflank/rflank = bp of the covered (M/=/X) interval abutting the gap on that side (0 at a contig
    end).  --rescueCompleteOnly uses these to demand a refmap alignment that reaches into BOTH flanks:
    a fill contiguous across gap+flanks HEALS the gap (recovers the sequence with no new path break),
    where an island fill -- covering only the gap middle -- leaves residual unrepresented query on each
    side and so adds a path fragment.  Same coverage definition as find_query_gaps."""
    covered = defaultdict(list)
    qlen = {}
    with _open(paf_path) as f:
        for line in f:
            if not line.strip() or line[0] == '@':
                continue
            t = line.rstrip("\n").split("\t")
            if len(t) < 12:
                continue
            q, ql, qs, qe, strand = t[0], int(t[1]), int(t[2]), int(t[3]), t[4]
            qlen[q] = ql
            cg = next((tag[5:] for tag in t[12:] if tag.startswith("cg:Z:")), None)
            if cg:
                covered[q].extend(_covered_query_intervals(qs, qe, strand, cg))
            else:
                covered[q].append((qs, qe))
    out = defaultdict(list)
    for q, ivs in covered.items():
        ivs.sort()
        merged = []
        for s, e in ivs:                           # union the aligned intervals
            if merged and s <= merged[-1][1]:
                merged[-1][1] = max(merged[-1][1], e)
            else:
                merged.append([s, e])
        pos, lflank = 0, 0                          # lflank = length of the covered run ending at pos
        for s, e in merged:
            if s - pos >= min_gap:
                out[q].append((pos, s, lflank, e - s))
            pos, lflank = e, e - s
        if qlen[q] - pos >= min_gap:               # trailing gap: no right flank (contig end)
            out[q].append((pos, qlen[q], lflank, 0))
    return out


def _spans_contiguous(lines, a, b):
    """True if the query spans of these refmap lines, unioned, cover [a,b] in a SINGLE gapless run --
    the reference alignment is contiguous across [a,b] (no hole), so a fill clipped from it reaches both
    flanks and heals the gap instead of islanding inside it."""
    ivs = sorted((int(l.split("\t")[2]), int(l.split("\t")[3])) for l in lines)
    cs = ce = None
    for s, e in ivs:
        if ce is None:
            cs, ce = s, e
        elif s <= ce:                              # overlap or abut -> same run
            ce = max(ce, e)
        else:
            if cs <= a and ce >= b:
                return True
            cs, ce = s, e
    return ce is not None and cs <= a and ce >= b


def find_insertion_segments(paf_path, backbone_nodes, min_seg):
    """{query: [(qs, qe, key)]} -- the query span of every minigraph record whose target node is NOT a
    reference-backbone node (a haplotype segment parked on an insertion node, a candidate detour),
    length >= min_seg.  `key` = (query, qs, qe, node) identifies the exact record so competitive
    re-anchoring can drop it from the pass-through if it re-anchors to the backbone.

    This is the OTHER half of the rescue's signal: find_query_gaps finds where the graph aligned
    NOTHING; this finds where it aligned to an insertion node instead of the reference it matches --
    the iterative-construction artifact (a genome that built its own node during construct and never
    aligned back to the identical backbone).  Backbone nodes are exactly the nodes the reference itself
    traverses (build_ref_node_table); anything else is rank>0 insertion sequence."""
    segs = defaultdict(list)
    with _open(paf_path) as f:
        for line in f:
            if not line.strip() or line[0] == '@':
                continue
            t = line.rstrip("\n").split("\t")
            if len(t) < 12 or t[5] in backbone_nodes:
                continue                               # already on the reference backbone
            qs, qe = int(t[2]), int(t[3])
            if qe - qs >= min_seg:
                segs[t[0]].append((qs, qe, (t[0], qs, qe, t[5])))
    return segs


def _record_identity(t):
    """Gap-excluded identity of an assembly->ref PAF record: 1 - de:f: if present, else the reported
    num_matches / block_length."""
    for tag in t[12:]:
        if tag.startswith("de:f:"):
            return max(0.0, min(1.0, 1.0 - float(tag[5:])))
    return int(t[9]) / (int(t[10]) or 1) if len(t) > 10 else 1.0


def _lookup_pid(pid_by_q, q, qs, qe):
    """Identity of the selected parent whose query interval contains this shattered piece's midpoint
    (1.0 if none -- shouldn't happen, since pieces come from shattering the selected parents)."""
    mid = (qs + qe) // 2
    for a, b, p in pid_by_q.get(q, []):
        if a <= mid < b:
            return p
    return 1.0


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
    # identity carried in via matches/block: gap_fill rewrites each shattered piece's num_matches to
    # its parent minimap2 record's real identity (paffy shatter itself sets num_matches == length, ie
    # 100%).  preserving that ratio through the node split makes every fill report honest identity, so
    # minIdentity can actually tell a good fill from a bad one.
    pid = int(t[9]) / (int(t[10]) or 1) if len(t) > 10 else 1.0
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
        matches = max(1, round(seg * pid))
        # AS: length-proportional score proxy (shatter drops the real minimap2 AS) that keeps the fill
        # above minScore so filter_paf judges on merit, not a defaulted AS:i:0.  num_matches + de:f:
        # carry the real identity so minIdentity applies.  rs:i:1 marks the record as a rescue fill.
        out.append("\t".join([q, qlen, str(q_s), str(q_e), strand, node, str(nlen),
                              str(n_s), str(n_e), str(matches), str(seg), mapq,
                              "tp:A:P", "AS:i:{}".format(seg), "rs:i:1",
                              "de:f:{:.4g}".format(1.0 - pid), "cg:Z:{}M".format(seg)]) + "\n")
    return out


def clip_piece_to_gaps(t, gaps):
    """Clip a gapless node-targeted PAF record (from lift_gapless_piece) to the query's gap
    intervals, so fills never overhang into minigraph-covered flanks (an overhang would give the
    contig both a forward and a reverse mapping to the same node -- a palindrome rgfa-split axes).
    Linear 1:1; handles reverse (query up => node down)."""
    qs, qe, strand = int(t[2]), int(t[3]), t[4]
    ns, ne = int(t[7]), int(t[8])
    pid = int(t[9]) / (int(t[10]) or 1)          # identity carried in from the lifted piece
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
        matches = max(1, round(seg * pid))
        # this clipped record is the final fill written to the merged PAF, so the honest score (AS),
        # identity (num_matches + de:f:) and rescue marker (rs:i:1) all have to live here too
        out.append("\t".join([t[0], t[1], str(iq_s), str(iq_e), strand, t[5], t[6],
                              str(n_s), str(n_e), str(matches), str(seg), t[11],
                              "tp:A:P", "AS:i:{}".format(seg), "rs:i:1",
                              "de:f:{:.4g}".format(1.0 - pid), "cg:Z:{}M".format(seg)]) + "\n")
    return out


# --- dipcall-style both-haplotype single-coverage confidence filter -------------------------------
# dipcall (github.com/lh3/dipcall) stays clean in SD by trusting a locus only where BOTH haplotypes
# align uniquely: `samflt` keeps primary, mapQ>=5, long alignments, and the confident BED is the
# single-coverage `paftools call` R regions INTERSECTED across haplotypes (`bedtk isec`).  We rebuild
# that here from the rescue's own asm->ref refmap and drop fills whose reference locus is not both-hap
# unique -- the SD-paralog collateral.  Validated on the 3 inversion carriers: the reconstruction
# matches dipcall's BED at jaccard~0.98 and culls 57-65% of collateral for ~2% recovery loss, with the
# flagship inversion 100% retained (carrier-truth/dipfilter/FINDING.md).

def _event_of(query):
    """`id=<sample>.<hap>|<contig>` -> (sample, event); event = <sample>.<hap>.  Falls back to the
    whole id if it doesn't parse, so an odd name just forms its own (unpaired) single-event group."""
    ev = query.split("|", 1)[0]
    if ev.startswith("id="):
        ev = ev[3:]
    sample = ev.rsplit(".", 1)[0] if "." in ev else ev
    return sample, ev


def _merge_intervals(ivs):
    if not ivs:
        return []
    ivs = sorted(ivs)
    out = [list(ivs[0])]
    for s, e in ivs[1:]:
        if s <= out[-1][1]:
            out[-1][1] = max(out[-1][1], e)
        else:
            out.append([s, e])
    return [(s, e) for s, e in out]


def _single_coverage(intervals):
    """Sub-intervals covered by EXACTLY one of the [s,e) intervals -- paftools `call` single-coverage
    `R` regions, via a coverage sweep line (emit each maximal run where the running depth is 1)."""
    pts = []
    for s, e in intervals:
        if e > s:
            pts.append((s, 1))
            pts.append((e, -1))
    pts.sort()
    out, cov, prev, i, n = [], 0, None, 0, len(pts)
    while i < n:
        pos = pts[i][0]
        if prev is not None and pos > prev and cov == 1:
            out.append((prev, pos))
        while i < n and pts[i][0] == pos:          # apply every delta at this coordinate together
            cov += pts[i][1]
            i += 1
        prev = pos
    return _merge_intervals(out)


def _intersect_intervals(a, b):
    """Intersect two sorted, non-overlapping interval lists."""
    out, i, j = [], 0, 0
    while i < len(a) and j < len(b):
        s, e = max(a[i][0], b[j][0]), min(a[i][1], b[j][1])
        if e > s:
            out.append((s, e))
        if a[i][1] <= b[j][1]:
            i += 1
        else:
            j += 1
    return out


def _overlap_len(starts, ivs, s, e):
    """Total overlap of [s,e) with sorted non-overlapping ivs (starts == [iv[0] for iv in ivs])."""
    if e <= s or not ivs:
        return 0
    i = max(0, bisect.bisect_right(starts, s) - 1)
    ov = 0
    while i < len(ivs) and ivs[i][0] < e:
        lo, hi = max(ivs[i][0], s), min(ivs[i][1], e)
        if hi > lo:
            ov += hi - lo
        i += 1
    return ov


def build_confident_intervals(rm_by_q, min_mapq=5, min_block=0, min_frac_samples=None):
    """Per sample, the both-haplotype single-coverage intersection of primary refmap alignments in
    reference coords -- dipcall's confident set, rebuilt from the rescue's own asm->ref refmap.
    Returns {sample: {ref_target: (starts, [(s,e)])}} for fast _overlap_len lookup.  A locus survives
    for a sample only if EVERY haplotype/event of that sample maps it with exactly one primary
    (mapQ>=min_mapq, block>=min_block) alignment, so SD paralogs (multi-coverage or low-mapQ in some
    hap) drop out.  If min_frac_samples is set (0..1), the per-sample sets are replaced by a single
    cohort mask kept where >= that fraction of samples are both-hap unique -- the population
    generalization of dipcall's diploid rule for non-diploid / multi-sample input; the same mask is
    then returned for every sample so the piece filter is unchanged."""
    by_event = defaultdict(lambda: defaultdict(list))     # event -> ref_target -> [(ts,te)]
    event_sample = {}
    for q, recs in rm_by_q.items():
        sample, event = _event_of(q)
        event_sample[event] = sample
        for qs, qe, tp, target, line in recs:
            if tp not in ('P', 'I'):
                continue
            t = line.rstrip("\n").split("\t")
            if len(t) < 12 or int(t[11]) < min_mapq or int(t[10]) < min_block:
                continue
            by_event[event][target].append((int(t[7]), int(t[8])))
    event_single = {ev: {c: _single_coverage(ivs) for c, ivs in chroms.items()}
                    for ev, chroms in by_event.items()}
    sample_events = defaultdict(list)
    for ev, sample in event_sample.items():
        sample_events[sample].append(ev)
    # per-sample both-hap (all-event) intersection
    sample_conf = {}                                      # sample -> {chrom: [(s,e)]}
    for sample, events in sample_events.items():
        cc = {}
        chroms = set()
        for ev in events:
            chroms |= set(event_single.get(ev, {}).keys())
        for c in chroms:
            acc = None
            for ev in events:
                s = event_single.get(ev, {}).get(c, [])
                if not s:                                # a hap with no unique mapping here
                    acc = None
                    break
                acc = s if acc is None else _intersect_intervals(acc, s)
                if not acc:
                    break
            if acc:
                cc[c] = acc
        sample_conf[sample] = cc
    if min_frac_samples is not None and sample_conf:
        cohort = _cohort_mask(sample_conf, min_frac_samples)
        sample_conf = {s: cohort for s in sample_conf}
    return {s: {c: ([iv[0] for iv in ivs], ivs) for c, ivs in cc.items()}
            for s, cc in sample_conf.items()}


def _cohort_mask(sample_conf, min_frac):
    """Cohort generalization: regions where >= min_frac of the samples are both-hap unique.  Per
    chrom, sweep every sample's confident intervals and keep the depth >= ceil(min_frac*N) runs."""
    n = len(sample_conf)
    need = max(1, (int(min_frac * n * 1000) + 999) // 1000)     # ceil(min_frac*n), no float/import
    by_chrom = defaultdict(list)
    for cc in sample_conf.values():
        for c, ivs in cc.items():
            by_chrom[c].extend(ivs)
    mask = {}
    for c, ivs in by_chrom.items():
        pts = []
        for s, e in ivs:
            if e > s:
                pts.append((s, 1)); pts.append((e, -1))
        pts.sort()
        out, cov, prev, i, m = [], 0, None, 0, len(pts)
        while i < m:
            pos = pts[i][0]
            if prev is not None and pos > prev and cov >= need:
                out.append((prev, pos))
            while i < m and pts[i][0] == pos:
                cov += pts[i][1]
                i += 1
            prev = pos
        merged = _merge_intervals(out)
        if merged:
            mask[c] = merged
    return mask


def _piece_confident(pt, confident, frac):
    """True if >= frac of a shattered asm->ref piece's reference span lies in its sample's both-hap
    confident set (pt is a PAF field list: pt[0]=query, pt[5]=ref_target, pt[7:9]=ref start/end)."""
    sample, _ = _event_of(pt[0])
    entry = confident.get(sample, {}).get(pt[5])
    if not entry:
        return False
    starts, ivs = entry
    ts, te = int(pt[7]), int(pt[8])
    return te > ts and _overlap_len(starts, ivs, ts, te) >= frac * (te - ts)


def gap_fill(minigraph_paf, rm_by_q, ref_table, out_paf, min_gap=1000, cover_frac=0.5,
             secondary_frac=0.5, paffy="paffy", assembly_fa=None, min_fill=0,
             confident_filter=False, confident_min_mapq=5, confident_min_block=0,
             confident_frac=0.5, confident_frac_samples=None, min_mapq=0, min_block=0,
             competitive=False, comp_min_mapq=30, comp_min_id=0.95, comp_min_block=0,
             complete_only=False, complete_flank=1000):
    """Write out_paf = the minigraph PAF (records pass through) + node-targeted reference records that
    fill its query coverage gaps.  rm_by_q is the ORIGINAL assembly->reference refmap (read_paf);
    ref_table is the reference->node table (build_ref_node_table).  min_mapq/min_block apply dipcall's
    samflt gate per source alignment in rescue_records; with confident_filter, fills whose reference
    locus is not in the both-haplotype single-coverage set are also dropped.

    With competitive=True, ALSO re-anchor: any graphmap record that parked a haplotype segment on an
    insertion node is DROPPED and replaced by a backbone anchor when the refmap maps that segment
    cleanly to the reference (rescue_records' uniqueness gate + the confident both-hap single-coverage
    copy-number guard + comp_min_mapq/comp_min_id).  This collapses the iterative-construction
    under-alignment (a genome stranded on its own redundant node); the confident filter is what keeps
    it from collapsing true insertions (a true insertion's ref locus is multi-covered -> not confident,
    so it survives).

    With complete_only=True, the gap path fills ONLY gaps a refmap alignment crosses contiguously into
    both flanks (find_flanked_gaps + _spans_contiguous): the fill abuts the flanking minigraph nodes and
    heals the gap, so recovery costs no new path fragment.  Island gaps -- where even the refmap breaks
    at the flanks -- are skipped.  Flank-anchoring is the confidence signal, so these fills bypass the
    min_block gate (which drops the very records that span the flanks), the confident filter, and the
    min_fill noise cut (their short node-boundary tiles ARE the heal; dropping them re-opens the gap).
    complete_flank = bp the alignment must reach into each flank.  Returns a stats dict."""
    gaps = None if complete_only else find_query_gaps(minigraph_paf, min_gap)
    confident = (build_confident_intervals(rm_by_q, confident_min_mapq, confident_min_block,
                                           confident_frac_samples)
                 if (confident_filter or competitive) else None)
    st = dict(gaps=0, filled=0, fill_records=0, unfilled=0, fill_filtered=0, detours=0, reanchored=0)
    # gate on the original refmap (tp:A: tags intact) -> the reference lines that cleanly fill gaps.
    # gap_selected is a flat pool of accepted primaries, so clip fills to the ACCEPTED gaps only: a
    # primary chosen for a clean gap can also span a REJECTED (ambiguous) gap of the same contig, and
    # clipping to all gaps would then fill the rejected one too -- reintroducing the SD collateral the
    # uniqueness gate just excluded.  (Safe: a primary overlapping another accepted gap in-gap would
    # have caused that gap's rejection, so accepted-only clipping never drops a legitimate fill.)
    gap_selected, accepted = [], defaultdict(list)
    if complete_only:
        # complete-only: fill a gap ONLY when a refmap alignment is contiguous across it AND reaches
        # >= complete_flank into BOTH minigraph-aligned flanks -- a fill that heals the gap instead of
        # islanding inside it.  Flank-anchoring is the confidence, so DROP the min_block gate here (it
        # kills the very records that span the flanks -- a ~31kb complete gap aligns a block under a
        # 50kb --rescueMinAlnBlock); rescue_records' uniqueness/secondary gates still apply.
        for q, ivs in find_flanked_gaps(minigraph_paf, min_gap).items():
            recs = rm_by_q.get(q, [])
            for gs, ge, lflank, rflank in ivs:
                st["gaps"] += 1
                fl, fr = min(complete_flank, lflank), min(complete_flank, rflank)
                resc = rescue_records(recs, gs, ge, cover_frac, secondary_frac, min_gap, min_mapq, 0)
                if resc and _spans_contiguous(resc, gs - fl, ge + fr):
                    gap_selected.extend(resc)
                    accepted[q].append((gs, ge))
                    st["filled"] += 1
                else:
                    st["unfilled"] += 1
    else:
        for q, ivs in gaps.items():
            for gs, ge in ivs:
                st["gaps"] += 1
                resc = rescue_records(rm_by_q.get(q, []), gs, ge, cover_frac, secondary_frac, min_gap,
                                      min_mapq, min_block)
                if resc:
                    gap_selected.extend(resc)
                    accepted[q].append((gs, ge))
                    st["filled"] += 1
                else:
                    st["unfilled"] += 1
    # competitive re-anchoring: the same gate, on insertion-node detours instead of gaps.  A detour whose
    # refmap alignment maps cleanly (unique, comp_min_mapq/id) to the backbone is only a CANDIDATE here;
    # whether its graphmap line is dropped is decided BELOW, once a surviving backbone fill actually
    # covers it -- so a candidate the confident copy-number guard or the node-boundary lift later kills
    # keeps its detour (never unanchored).  Backbone nodes = the nodes the reference itself traverses.
    comp_selected, candidates, detour_by_q = [], [], defaultdict(list)
    if competitive:
        backbone_nodes = {node for entries in ref_table.values()
                          for (_rqs, _rqe, node, _nl, _nts, _nte) in entries}
        for q, seglist in find_insertion_segments(minigraph_paf, backbone_nodes, min_gap).items():
            for qs, qe, key in seglist:
                st["detours"] += 1
                resc = rescue_records(rm_by_q.get(q, []), qs, qe, cover_frac, secondary_frac, min_gap,
                                      comp_min_mapq, comp_min_block)
                if not resc or any(_record_identity(l.rstrip("\n").split("\t")) < comp_min_id
                                   for l in resc):
                    continue                           # no clean unique backbone home for this detour
                comp_selected.extend(resc)
                candidates.append((key, q, qs, qe))
                detour_by_q[q].append((qs, qe))

    def _shatter_lift(selected, apply_confident, clip_regions, tag):
        """selected asm->ref lines -> node-targeted, query-clipped fills.  paffy shatter -> gapless
        pieces; each is (optionally) confident-filtered, its identity corrected back to its parent's
        (shatter emits every piece at 100%), lifted to node coords, and clipped to clip_regions.  Dedupe
        first: a parent spanning several targets is appended once each.  (fills, n_confident_filtered)."""
        out_f, nfilt = [], 0
        if not selected:
            return out_f, nfilt
        selected = list(dict.fromkeys(selected))
        sel_path = out_paf + "." + tag + ".paf"
        with open(sel_path, "w") as fh:
            fh.writelines(selected)
        pid_by_q, seen = defaultdict(list), set()
        for pline in selected:
            ppt = pline.split("\t")
            k = (ppt[0], ppt[2], ppt[3])
            if k not in seen:
                seen.add(k)
                pid_by_q[ppt[0]].append((int(ppt[2]), int(ppt[3]), _record_identity(ppt)))
        shattered = subprocess.run([paffy, "shatter", "-i", sel_path],
                                   capture_output=True, text=True, check=True).stdout
        for line in shattered.splitlines():
            pt = line.split("\t")
            if apply_confident and confident is not None and not _piece_confident(pt, confident, confident_frac):
                nfilt += 1                             # ref locus not both-hap unique -> SD-paralog collateral
                continue
            pt[9] = str(max(1, round(int(pt[10]) * _lookup_pid(pid_by_q, pt[0], int(pt[2]), int(pt[3])))))
            for piece in lift_gapless_piece(pt, ref_table):
                ppt = piece.rstrip("\n").split("\t")
                out_f.extend(clip_piece_to_gaps(ppt, clip_regions.get(ppt[0], [])))
        return out_f, nfilt

    # gap fills: confident-filter ONLY if --rescueConfidentFilter -- competitive alone must not change
    # the gap path (build_confident_intervals is built for competitive, but the gap fills don't consult it).
    # complete-only bypasses it too: flank-anchoring is the gate, and dropping interior pieces would
    # re-open a healed gap.
    gap_fills, nf = _shatter_lift(gap_selected, confident_filter and not complete_only, accepted, "sel")
    st["fill_filtered"] += nf
    # bounds safety net: never emit a fill whose node interval is out of bounds or whose query/node
    # spans disagree -- cactus_consolidated rejects such records and the whole run dies hours later.
    def _bounds_ok(c):
        nlen, n_s, n_e = int(c[6]), int(c[7]), int(c[8])
        return 0 <= n_s < n_e <= nlen and (int(c[3]) - int(c[2])) == (n_e - n_s)

    # gap fills: bounds net + the min_fill noise filter (drops tiny node-boundary fragments).
    gap_valid, dropped, short = [], 0, 0
    for line in gap_fills:
        c = line.rstrip("\n").split("\t")
        if not _bounds_ok(c):
            dropped += 1
        elif not complete_only and int(c[8]) - int(c[7]) < min_fill:   # tiny node-boundary fragment ->
            short += 1                                                 # drop as noise (complete-only's
        else:                                                          # tiles heal the gap: keep them)
            gap_valid.append(line)

    # competitive fills: copy-number guard ALWAYS on.  Remove a detour ONLY where a piece that will
    # actually be EMITTED covers >= cover_frac of it, and emit only pieces landing in a removed detour
    # (else the kept detour + fill would double-anchor the span).  Competitive pieces are deliberately
    # NOT min_fill-filtered: a re-anchor legitimately tiles many short backbone nodes, so applying the
    # gap path's per-piece min_fill would drop them AND -- because the detour is then removed on their
    # coverage -- strand the detour's whole span unanchored, which the clip phase deletes.  So gate the
    # removal and the emission on the same bounds-valid (min_fill-exempt) pieces.
    remove, comp_fills = set(), []
    if comp_selected:
        cf, nf = _shatter_lift(comp_selected, True, detour_by_q, "comp")
        st["fill_filtered"] += nf
        surviving = []
        for f in cf:
            if _bounds_ok(f.rstrip("\n").split("\t")):
                surviving.append(f)
            else:
                dropped += 1
        cov = defaultdict(list)
        for f in surviving:
            c = f.split("\t"); cov[c[0]].append((int(c[2]), int(c[3])))
        kept = defaultdict(list)
        for key, q, qs, qe in candidates:
            if covered_fraction(cov.get(q, []), qs, qe) >= cover_frac:
                remove.add(key); kept[q].append((qs, qe)); st["reanchored"] += 1
        for f in surviving:
            c = f.split("\t"); mid = (int(c[2]) + int(c[3])) // 2
            if any(a <= mid < b for a, b in kept.get(c[0], [])):
                comp_fills.append(f)
    fills = gap_valid + comp_fills
    st["fill_records"] = len(fills)
    st["fill_dropped"] = dropped
    st["fill_short"] = short
    with open(out_paf, "w") as out:
        with _open(minigraph_paf) as f:               # minigraph records pass through (minus re-anchored detours)
            for line in f:
                if remove:
                    t = line.rstrip("\n").split("\t")
                    if (len(t) >= 12 and t[2].isdigit() and t[3].isdigit()
                            and (t[0], int(t[2]), int(t[3]), t[5]) in remove):
                        continue                        # this detour was competitively re-anchored to the backbone
                out.write(line if line.endswith("\n") else line + "\n")
        for line in fills:
            out.write(line)
    return st


def gap_fill_pafs(minigraph_paf, refmap_paf, out_paf, ref_prefix, min_gap=1000, cover_frac=0.5,
                  secondary_frac=0.5, paffy="paffy", assembly_fa=None, min_fill=0,
                  confident_filter=False, confident_min_mapq=5, confident_min_block=0,
                  confident_frac=0.5, confident_frac_samples=None, min_mapq=0, min_block=0,
                  competitive=False, comp_min_mapq=30, comp_min_id=0.95, comp_min_block=0,
                  complete_only=False, complete_flank=1000):
    """Convenience: build the ref->node table from the minigraph PAF and gap-fill one file."""
    table = build_ref_node_table(minigraph_paf, ref_prefix)
    return gap_fill(minigraph_paf, read_paf(refmap_paf), table, out_paf, min_gap, cover_frac,
                    secondary_frac, paffy, assembly_fa, min_fill, confident_filter,
                    confident_min_mapq, confident_min_block, confident_frac, confident_frac_samples,
                    min_mapq, min_block, competitive, comp_min_mapq, comp_min_id, comp_min_block,
                    complete_only, complete_flank)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--graphmap", required=True, help="minigraph graphmap PAF (may be .gz)")
    ap.add_argument("--refmap", required=True, help="assembly->reference PAF from cactus-refmap (may be .gz)")
    ap.add_argument("-o", "--output", required=True, help="output merged PAF")
    ap.add_argument("--ref-prefix", required=True, help="reference query prefix, e.g. 'id=CHM13|'")
    ap.add_argument("--min-gap", type=int, default=1000)
    ap.add_argument("--cover-frac", type=float, default=0.5)
    ap.add_argument("--secondary-frac", type=float, default=0.5)
    ap.add_argument("--min-fill", type=int, default=0, help="drop fills shorter than this many node bp")
    ap.add_argument("--paffy", default="paffy", help="path to the paffy binary [paffy]")
    ap.add_argument("--assembly", default=None, help="query fasta, to also surface fully-missing contigs (paffy -q)")
    # dipcall-style both-haplotype single-coverage confidence filter (SD-paralog collateral cull)
    ap.add_argument("--confident-filter", action="store_true",
                    help="drop fills whose reference locus is not both-haplotype single-coverage unique")
    ap.add_argument("--confident-min-mapq", type=int, default=5,
                    help="min primary mapQ for an alignment to count toward the confident set [5]")
    ap.add_argument("--confident-min-block", type=int, default=0,
                    help="min alignment block length to count toward the confident set [0]")
    ap.add_argument("--confident-frac", type=float, default=0.5,
                    help="min fraction of a fill piece that must lie in the confident set to keep it [0.5]")
    ap.add_argument("--confident-frac-samples", type=float, default=None,
                    help="cohort generalization: keep loci where >= this fraction of samples are "
                         "both-hap unique, instead of per-sample both-hap (for non-diploid input)")
    # dipcall samflt gate: only fill a gap from source alignments >= these mapQ / block-length
    ap.add_argument("--min-mapq", type=int, default=0,
                    help="only fill from source refmap alignments with mapQ >= this (dipcall samflt: 5)")
    ap.add_argument("--min-aln-block", type=int, default=0,
                    help="only fill from source refmap alignments with >= this bp of block (dipcall: 50000)")
    # competitive re-anchoring: move insertion-node detours onto the reference backbone
    ap.add_argument("--competitive", action="store_true",
                    help="re-anchor insertion-node detours to the reference where the refmap maps them "
                         "cleanly (collapses iterative-construction under-alignment; confident guard forced on)")
    ap.add_argument("--competitive-min-mapq", type=int, default=30,
                    help="min source mapQ for a competitive backbone re-anchor [30]")
    ap.add_argument("--competitive-min-id", type=float, default=0.95,
                    help="min source identity for a competitive backbone re-anchor [0.95]")
    ap.add_argument("--competitive-min-block", type=int, default=0,
                    help="min source block length for a competitive backbone re-anchor [0]")
    # complete-only: fill only gaps a refmap alignment crosses contiguously into both flanks
    ap.add_argument("--complete-only", action="store_true",
                    help="fill only gaps a refmap alignment crosses contiguously into both flanks "
                         "(heals the gap, no new path fragment); skip island gaps.  min_block / min_fill "
                         "exempt, confident filter bypassed -- flank-anchoring is the gate")
    ap.add_argument("--complete-flank", type=int, default=1000,
                    help="--complete-only: bp a refmap alignment must reach into each flank [1000]")
    args = ap.parse_args()
    st = gap_fill_pafs(args.graphmap, args.refmap, args.output, args.ref_prefix, args.min_gap,
                       args.cover_frac, args.secondary_frac, args.paffy, args.assembly, args.min_fill,
                       args.confident_filter, args.confident_min_mapq, args.confident_min_block,
                       args.confident_frac, args.confident_frac_samples, args.min_mapq, args.min_aln_block,
                       args.competitive, args.competitive_min_mapq, args.competitive_min_id,
                       args.competitive_min_block, args.complete_only, args.complete_flank)
    sys.stderr.write("[rescue] {gaps} gaps: {filled} filled -> {fill_records} node records, "
                     "{unfilled} unfilled, {fill_filtered} pieces dropped by confidence filter; "
                     "{detours} detours: {reanchored} re-anchored to backbone\n".format(**st))


if __name__ == "__main__":
    main()
