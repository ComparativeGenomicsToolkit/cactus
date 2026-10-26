#!/usr/bin/env python3
"""Score a cactus alignment against BUSCO gene models.

For each ordered pair of genomes (A, B), take the BUSCO genes that compleasm calls Single in both,
with no frameshift and a CDS a whole number of codons long in both: a 1:1 ortholog with known
exon structure on each side.  Every base of each such gene in A, plus FLANK bp either side, is
lifted over to B through the HAL (halLiftover --outPSL) and classed by where it lands:

  coding bases (A's CDS):
    in phase      in B's CDS of the same gene, at the same codon position.  Correct alignment of
                  a coding base does this; misplacing it within the right exon does so only 1 time
                  in 3, and landing anywhere else never does
    out of phase  in B's CDS of the same gene, at a different codon position
    in locus      in B's copy of the gene (+-MARGIN bp) but not its CDS
    elsewhere     aligned, but outside B's copy of the gene -- certainly wrong
    unaligned     not aligned to B at all
  intron and flank bases (A's gene outside its CDS, and FLANK bp either side):
    in locus / elsewhere / unaligned, as above.  There is no base-level truth outside the CDS, but
    anything that lands outside the ortholog's locus is certainly wrong, which is where
    alignment made from noise shows up.

A base aligned to several places in B (a duplication) is classed by the best of them, and the
share of aligned bases with more than one landing is reported.

Genes come from compleasm's full_table.tsv (found as BUSCO_DIR/<genome>/**/full_table.tsv): the
Codons column gives CDS exons, 1-based inclusive, in transcript order.  Sequence names must match
the HAL's; the script stops if fewer than 90% of a genome's BUSCO sequences are in the HAL.

usage: busco-align-eval.py --hal ALN.hal --busco BUSCO_DIR [--flank 5000] [--margin 20000]
                           [--threads 8] GENOME1 GENOME2 ...
       busco-align-eval.py --summarize [NAME1=]arm1/eval.txt [NAME2=]arm2/eval.txt ...
"""
import sys, os, glob, argparse, subprocess, tempfile, collections
from multiprocessing import Pool
import numpy as np

CDS_CLASSES = ['in phase', 'out of phase', 'in locus', 'elsewhere', 'unaligned']
NC_CLASSES = ['in locus', 'elsewhere', 'unaligned']
# best-of ranking per base: higher wins
R_UNALIGNED, R_ELSEWHERE, R_LOCUS, R_OUTPHASE, R_INPHASE = 0, 1, 2, 3, 4

Gene = collections.namedtuple('Gene', 'seq start end strand cds')   # 0-based half-open; cds [(s, e)] in transcript order

def read_busco(path):
    """BUSCO id -> Gene, for Single genes with no frameshift and a whole number of codons"""
    genes = {}
    with open(path) as f:
        header = f.readline().rstrip('\n').split('\t')
        col = {h: i for i, h in enumerate(header)}
        for line in f:
            p = line.rstrip('\n').split('\t')
            if p[col['Status']] != 'Single' or p[col['Frameshift events']] != '0':
                continue
            cds = []
            for b in p[col['Codons']].split('|'):
                s, e, st = b.split('_')
                cds.append((int(s) - 1, int(e)))
            if sum(e - s for s, e in cds) % 3:
                continue
            seq = p[col['Sequence']]
            genes[p[col['Gene']]] = Gene(seq, min(s for s, e in cds), max(e for s, e in cds), p[col['Strand']], cds)
    return genes

def find_table(busco_dir, genome):
    hits = glob.glob(os.path.join(busco_dir, genome, '**', 'full_table.tsv'), recursive=True)
    if len(hits) != 1:
        sys.exit('expected one full_table.tsv under {}, found {}'.format(os.path.join(busco_dir, genome), hits))
    return hits[0]

def chrom_sizes(hal, genome):
    out = subprocess.run(['halStats', '--chromSizes', genome, hal], check=True, capture_output=True, text=True).stdout
    return {l.split('\t')[0]: int(l.split('\t')[1]) for l in out.splitlines() if l.strip()}

def codon_array(g, lo, hi):
    """codon position (1-3) of each base in [lo, hi) that is in g's CDS, else 0"""
    a = np.zeros(hi - lo, dtype=np.int8)
    k = 0
    for s, e in g.cds:
        n = e - s
        pos = (np.arange(k, k + n) % 3 + 1).astype(np.int8)
        if g.strand == '+':
            a[s - lo:e - lo] = pos
        else:
            a[s - lo:e - lo] = pos[::-1]      # transcript runs from e-1 down to s
        k += n
    return a

def score_pair(args):
    hal, a_name, b_name, genes_a, genes_b, sizes_a, flank, margin = args
    shared = sorted(set(genes_a) & set(genes_b))
    # A regions: gene +- flank, clipped to the sequence
    regions = {}
    with tempfile.NamedTemporaryFile('w', suffix='.bed', delete=False) as bed:
        for gid in shared:
            g = genes_a[gid]
            if g.seq not in sizes_a:
                continue
            lo, hi = max(0, g.start - flank), min(sizes_a[g.seq], g.end + flank)
            regions[gid] = (lo, hi)
            bed.write('{}\t{}\t{}\t{}\n'.format(g.seq, lo, hi, gid))
        bed_path = bed.name
    psl_path = bed_path + '.psl'
    subprocess.run(['halLiftover', '--outPSLWithName', hal, a_name, bed_path, b_name, psl_path],
                   check=True, stderr=subprocess.DEVNULL)

    # per gene: A class of each region base (0 flank, 1 intron, 2 CDS) with codon position, and
    # the best landing and landing count seen so far
    state = {}
    for gid, (lo, hi) in regions.items():
        g = genes_a[gid]
        cpos = codon_array(g, lo, hi)
        kind = np.zeros(hi - lo, dtype=np.int8)
        kind[g.start - lo:g.end - lo] = 1
        kind[cpos > 0] = 2
        state[gid] = [lo, cpos, kind, np.zeros(hi - lo, dtype=np.int8), np.zeros(hi - lo, dtype=np.int32)]
    # B: the ortholog's locus and codon positions over it
    locus_b = {}
    for gid in regions:
        g = genes_b[gid]
        lo, hi = max(0, g.start - margin), g.end + margin
        locus_b[gid] = (g.seq, lo, hi, codon_array(g, lo, hi))

    with open(psl_path) as f:
        for line in f:
            p = line.rstrip('\n').split('\t')
            gid, p = p[0], p[1:]
            if gid not in state:
                continue
            strand, t_name, t_size = p[8], p[13], int(p[14])
            sizes = [int(x) for x in p[18].rstrip(',').split(',')]
            q_starts = [int(x) for x in p[19].rstrip(',').split(',')]
            t_starts = [int(x) for x in p[20].rstrip(',').split(',')]
            lo, cpos, kind, best, count = state[gid]
            b_seq, b_lo, b_hi, b_cpos = locus_b[gid]
            for n, qs, ts in zip(sizes, q_starts, t_starts):
                qi = np.arange(qs, qs + n) - lo
                keep = (qi >= 0) & (qi < len(best))
                if not keep.any():
                    continue
                if strand[-1] == '-':      # block target coordinates are on B's reverse strand
                    tpos = t_size - 1 - (ts + np.arange(n))
                else:
                    tpos = ts + np.arange(n)
                qi, tpos = qi[keep], tpos[keep]
                rank = np.full(len(qi), R_ELSEWHERE, dtype=np.int8)
                if t_name == b_seq:
                    inl = (tpos >= b_lo) & (tpos < b_hi)
                    rank[inl] = R_LOCUS
                    bc = np.zeros(len(qi), dtype=np.int8)
                    bc[inl] = b_cpos[tpos[inl] - b_lo]
                    coding = (kind[qi] == 2) & (bc > 0)
                    rank[coding] = np.where(bc[coding] == cpos[qi][coding], R_INPHASE, R_OUTPHASE)
                np.maximum.at(best, qi, rank)
                np.add.at(count, qi, 1)
    os.unlink(bed_path)
    os.unlink(psl_path)

    tally = {'cds': np.zeros(5, np.int64), 'intron': np.zeros(5, np.int64), 'flank': np.zeros(5, np.int64)}
    multi = {'cds': 0, 'intron': 0, 'flank': 0}
    for gid, (lo, cpos, kind, best, count) in state.items():
        for k, name in ((2, 'cds'), (1, 'intron'), (0, 'flank')):
            m = kind == k
            tally[name] += np.bincount(best[m], minlength=5)[:5]
            multi[name] += int((count[m] > 1).sum())
    return a_name, b_name, len(regions), tally, multi

def pct(x, tot):
    return 100.0 * x / tot if tot else 0.0

def report(results, out):
    out.write('# BUSCO-based alignment scoring.  Rows are ordered pairs A -> B; bases are A\'s.\n')
    out.write('# coding: share of CDS bases in each class; "precision" = in phase / aligned\n')
    out.write('{:36s} {:>6s} {:>12s} {:>8s} {:>8s} {:>8s} {:>8s} {:>8s} {:>9s} {:>6s}\n'.format(
        'pair', 'genes', 'CDS bp', 'inphase', 'outph', 'inlocus', 'elsewh', 'unalig', 'precision', 'multi'))
    for a, b, n, t, m in results:
        c = t['cds']; tot = c.sum(); aligned = tot - c[R_UNALIGNED]
        out.write('{:36s} {:>6} {:>12,} {:>7.1f}% {:>7.1f}% {:>7.1f}% {:>7.1f}% {:>7.1f}% {:>8.1f}% {:>5.1f}%\n'.format(
            a + ' -> ' + b, n, int(tot), pct(c[R_INPHASE], tot), pct(c[R_OUTPHASE], tot), pct(c[R_LOCUS], tot),
            pct(c[R_ELSEWHERE], tot), pct(c[R_UNALIGNED], tot), pct(c[R_INPHASE], aligned), pct(m['cds'], aligned)))
    for part in ('intron', 'flank'):
        out.write('\n# {}: share of bases in each class; "in-locus share" = in locus / aligned\n'.format(part))
        out.write('{:36s} {:>14s} {:>8s} {:>8s} {:>8s} {:>13s} {:>6s}\n'.format(
            'pair', 'bp', 'inlocus', 'elsewh', 'unalig', 'in-locus share', 'multi'))
        for a, b, n, t, m in results:
            c = t[part]; tot = c.sum(); aligned = tot - c[R_UNALIGNED]
            out.write('{:36s} {:>14,} {:>7.1f}% {:>7.1f}% {:>7.1f}% {:>12.1f}% {:>5.1f}%\n'.format(
                a + ' -> ' + b, int(tot), pct(c[R_LOCUS], tot), pct(c[R_ELSEWHERE], tot),
                pct(c[R_UNALIGNED], tot), pct(c[R_LOCUS], aligned), pct(m[part], aligned)))
    # machine-readable copy for --summarize
    for a, b, n, t, m in results:
        for part in ('cds', 'intron', 'flank'):
            out.write('#TSV\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(a, b, n, part, ','.join(str(int(x)) for x in t[part]), m[part]))

def summarize(paths, out):
    """one row per pair, one column group per arm"""
    data = collections.OrderedDict()
    arms = []
    for path in paths:
        # NAME=PATH names the arm NAME; a bare PATH is named by its directory
        arm, eq, path = path.rpartition('=')
        if not eq:
            arm = os.path.basename(os.path.dirname(os.path.abspath(path)))
        arms.append(arm)
        for line in open(path):
            if line.startswith('#TSV\t'):
                _, a, b, n, part, counts, multi = line.rstrip('\n').split('\t')
                data.setdefault((a, b), {}).setdefault(arm, {})[part] = np.array([int(x) for x in counts.split(',')])
    out.write('# aligned = share of A\'s bases aligned to B; precision = coding bases in phase / aligned;\n'
              '# in-locus = intron+flank bases landing in the ortholog locus / aligned\n')
    head = '{:36s}'.format('pair') + ''.join(' | {:^26s}'.format(a) for a in arms)
    out.write(head + '\n' + '{:36s}'.format('') + ''.join(' | {:>8s} {:>8s} {:>8s}'.format('CDS aln', 'precis', 'nc inloc') for a in arms) + '\n')
    for (a, b), per_arm in data.items():
        row = '{:36s}'.format(a + ' -> ' + b)
        for arm in arms:
            d = per_arm.get(arm)
            if d is None:
                row += ' | {:>26s}'.format('-')
                continue
            c = d['cds']; tot = c.sum(); al = tot - c[R_UNALIGNED]
            nc = d['intron'] + d['flank']; ntot = nc.sum(); nal = ntot - nc[R_UNALIGNED]
            row += ' | {:>7.1f}% {:>7.1f}% {:>7.1f}%'.format(pct(al, tot), pct(c[R_INPHASE], al), pct(nc[R_LOCUS], nal))
        out.write(row + '\n')

def main():
    if len(sys.argv) > 1 and sys.argv[1] == '--summarize':
        summarize(sys.argv[2:], sys.stdout)
        return
    ap = argparse.ArgumentParser()
    ap.add_argument('--hal', required=True)
    ap.add_argument('--busco', required=True)
    ap.add_argument('--flank', type=int, default=5000)
    ap.add_argument('--margin', type=int, default=20000)
    ap.add_argument('--threads', type=int, default=8)
    ap.add_argument('genomes', nargs='+')
    o = ap.parse_args()

    genes, sizes = {}, {}
    for g in o.genomes:
        genes[g] = read_busco(find_table(o.busco, g))
        sizes[g] = chrom_sizes(o.hal, g)
        seqs = {x.seq for x in genes[g].values()}
        found = len([s for s in seqs if s in sizes[g]])
        sys.stderr.write('{}: {} usable Single BUSCOs on {} sequences, {} of them in the HAL\n'.format(
            g, len(genes[g]), len(seqs), found))
        if seqs and found < 0.9 * len(seqs):
            sys.exit('{}: BUSCO sequence names do not match the HAL\'s (e.g. {} vs {})'.format(
                g, sorted(seqs)[:3], sorted(sizes[g])[:3]))
    jobs = [(o.hal, a, b, genes[a], genes[b], sizes[a], o.flank, o.margin)
            for a in o.genomes for b in o.genomes if a != b]
    with Pool(o.threads) as pool:
        results = pool.map(score_pair, jobs)
    report(results, sys.stdout)

if __name__ == '__main__':
    main()
