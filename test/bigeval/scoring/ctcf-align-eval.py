#!/usr/bin/env python3
"""Score a cactus alignment by whether it carries human CTCF sites to CTCF motifs in other genomes.

The idea is Ofir Raz's -- does an aligned regulatory element keep its transcription-factor motif
in the other species -- made into a per-alignment number with its own chance baseline.

Sites.  ENCODE cCREs of the chosen classes (default CA-CTCF: accessible chromatin with CTCF signal)
in the source genome (human, GRCh38 coordinates), each reduced to its single best strong CTCF
motif hit, if it has one.  The motif instance, not the whole cCRE, is what is lifted, so where it
lands can be scored at motif resolution.

Motif hits.  JASPAR MA0139.1, scanned genome-wide on both strands in every genome once (scan
subcommand).  "Strong" is an exact p-value against a uniform background, from the integer score
distribution of the PWM (as FIMO does): p <= 1e-4 by default, FIMO's own default.  Stricter is
cleaner but thin: only ~6% of V4 CA-CTCF cCREs have a p <= 1e-5 hit (~20% at 1e-4), because V4 calls
a cCRE CTCF on CTCF signal in any of hundreds of biosamples.

Scoring (score subcommand), per target genome:
  aligned     any base of the source motif is aligned to the target (halLiftover)
  conserved   a strong motif hit in the target at the landing -- overlapping the landing's span
              +-PAD by at least half the motif -- on the strand the alignment carries the source
              motif's strand to.  Judged at the primary landing (the one with most aligned bases)
              and, separately, at any landing.
  background  the same test, same strand, at 20 windows shifted 1-10 kb either side of the primary
              landing: the chance rate in that neighbourhood's own sequence.  expected = aligned x
              background rate, and excess = conserved (primary) - expected.

An alignment that carries real orthologous CTCF sites shows conserved well above background; one
that aligns sequence to non-homologous places shows conserved at background, however much it
aligns.  summarize compares arms, and in particular what each arm aligns that a baseline arm
does not -- Ofir's question: does the extra alignment carry the motif?

usage:
  ctcf-align-eval.py scan --fasta GENOME.fa[.gz] --out HITS.npz [--pvalue 1e-4]
  ctcf-align-eval.py score --hal ALN.hal --source GENOME --ccre cCREs.bed --hits DIR
                           [--classes CA-CTCF] [--pvalue 1e-4] [--threads 4] [--sites OUT.tsv.gz]
                           TARGET1 TARGET2 ...
  ctcf-align-eval.py summarize [--baseline ARM] [ARM1=]ARM1/ctcf-sites.tsv.gz [ARM2=]ARM2/ctcf-sites.tsv.gz ...

DIR holds one GENOME.npz per genome from scan, named by its HAL genome name.
"""
import sys, os, gzip, argparse, subprocess, tempfile, collections
from multiprocessing import Pool
import numpy as np

# JASPAR MA0139.1, CTCF: counts, rows A C G T, 19 columns
MA0139_1 = np.array([
    [87, 167, 281, 56, 8, 744, 40, 107, 851, 5, 333, 54, 12, 56, 104, 372, 82, 117, 402],
    [291, 145, 49, 800, 903, 13, 528, 433, 11, 0, 3, 12, 0, 8, 733, 13, 482, 322, 181],
    [76, 414, 449, 21, 0, 65, 334, 48, 32, 903, 566, 504, 890, 775, 5, 507, 307, 73, 266],
    [459, 187, 134, 36, 2, 91, 11, 324, 18, 3, 9, 341, 8, 71, 67, 17, 37, 396, 59]], dtype=float)
SCALE = 100            # log2-odds are kept as integers, x SCALE
N_PENALTY = -10 ** 6   # any N in a window rules it out
LUT = np.full(256, 4, dtype=np.uint8)
for _c, _v in zip(b'ACGTacgt', [0, 1, 2, 3, 0, 1, 2, 3]):
    LUT[_c] = _v
CONTROL_OFFSETS = [s * k * 1000 for k in range(1, 11) for s in (1, -1)]

def pwm():
    """forward and reverse-complement integer score matrices, each w x 5 (A C G T N)"""
    counts = MA0139_1
    freq = (counts + 0.25) / (counts.sum(axis=0) + 1.0)
    lod = np.round(np.log2(freq / 0.25) * SCALE).astype(np.int64)        # 4 x w
    w = lod.shape[1]
    fwd = np.full((w, 5), N_PENALTY, dtype=np.int64)
    fwd[:, :4] = lod.T
    rev = np.full((w, 5), N_PENALTY, dtype=np.int64)
    rev[:, :4] = lod.T[::-1, ::-1]       # reversed columns, complemented bases (A C G T -> T G C A)
    return fwd, rev

def threshold(fwd, pvalue):
    """smallest integer score whose upper tail has probability <= pvalue under a uniform background"""
    dist, offset = np.ones(1), 0
    for row in fwd[:, :4]:
        lo, hi = row.min(), row.max()
        new = np.zeros(len(dist) + hi - lo)
        for v in row:
            new[v - lo:v - lo + len(dist)] += 0.25 * dist
        dist, offset = new, offset + lo
    tail = np.cumsum(dist[::-1])[::-1]       # tail[k] = P(score >= offset + k)
    return int(offset + np.argmax(tail <= pvalue))

def opener(path):
    return gzip.open(path, 'rb') if path.endswith('.gz') else open(path, 'rb')

def read_fasta(path):
    name, chunks = None, []
    with opener(path) as f:
        for line in f:
            if line.startswith(b'>'):
                if name is not None:
                    yield name, b''.join(chunks)
                name, chunks = line[1:].split()[0].decode(), []
            else:
                chunks.append(line.rstrip())
    if name is not None:
        yield name, b''.join(chunks)

def scan_sequence(seq, fwd, rev, thr, chunk=20_000_000):
    """(starts, strands, scores) of all windows scoring >= thr, strands +1/-1"""
    codes = LUT[np.frombuffer(seq, dtype=np.uint8)]
    w = fwd.shape[0]
    out_s, out_d, out_v = [], [], []
    for cs in range(0, max(0, len(codes) - w + 1), chunk):
        n = min(chunk, len(codes) - w + 1 - cs)
        sf = np.zeros(n, dtype=np.int64)
        sr = np.zeros(n, dtype=np.int64)
        for j in range(w):
            c = codes[cs + j:cs + j + n]
            sf += fwd[j][c]
            sr += rev[j][c]
        for s, d in ((sf, 1), (sr, -1)):
            i = np.flatnonzero(s >= thr)
            out_s.append(i + cs)
            out_d.append(np.full(len(i), d, dtype=np.int8))
            out_v.append(s[i])
    if not out_s:
        return np.zeros(0, np.int64), np.zeros(0, np.int8), np.zeros(0, np.int64)
    s, d, v = np.concatenate(out_s), np.concatenate(out_d), np.concatenate(out_v)
    o = np.argsort(s, kind='stable')
    return s[o], d[o], v[o]

def cmd_scan(args):
    fwd, rev = pwm()
    thr = threshold(fwd, args.pvalue)
    names, seqi, starts, strands, scores, lengths = [], [], [], [], [], []
    for name, seq in read_fasta(args.fasta):
        s, d, v = scan_sequence(seq, fwd, rev, thr)
        seqi.append(np.full(len(s), len(names), dtype=np.int32))
        names.append(name)
        lengths.append(len(seq))
        starts.append(s); strands.append(d); scores.append(v)
    np.savez(args.out, names=np.array(names), lengths=np.array(lengths, dtype=np.int64),
             seq=np.concatenate(seqi), start=np.concatenate(starts), strand=np.concatenate(strands),
             score=np.concatenate(scores), threshold=thr, width=fwd.shape[0], pvalue=args.pvalue)
    total = sum(lengths)
    sys.stderr.write('{}: {:,} bp in {:,} sequences, {:,} hits at p <= {:g} (score >= {}); chance alone '
                     'would give ~{:,}\n'.format(args.fasta, total, len(names), sum(len(x) for x in starts),
                                                 args.pvalue, thr, int(2 * total * args.pvalue)))

class Hits:
    """a genome's motif hits, per sequence, sorted by start"""
    def __init__(self, path, pvalue):
        z = np.load(path)
        if abs(float(z['pvalue']) - pvalue) > 1e-12:
            sys.exit('{} was scanned at p <= {}, not {}'.format(path, float(z['pvalue']), pvalue))
        self.width = int(z['width'])
        self.lengths = dict(zip(z['names'].tolist(), z['lengths'].tolist()))
        self.by_seq = {}
        seq, start, strand, score = z['seq'], z['start'], z['strand'], z['score']
        order = np.lexsort((start, seq))
        seq, start, strand, score = seq[order], start[order], strand[order], score[order]
        bounds = np.flatnonzero(np.diff(seq)) + 1
        for lo, hi in zip(np.concatenate(([0], bounds)), np.concatenate((bounds, [len(seq)]))):
            if hi > lo:
                self.by_seq[z['names'][seq[lo]]] = (start[lo:hi], strand[lo:hi], score[lo:hi])

    def has_hit(self, name, lo, hi, strand, min_overlap):
        """is there a hit on strand overlapping [lo, hi) by at least min_overlap bases"""
        rec = self.by_seq.get(name)
        if rec is None:
            return False
        starts, strands, _ = rec
        w = self.width
        i = int(np.searchsorted(starts, lo - w + min_overlap))
        while i < len(starts) and starts[i] <= hi - min_overlap:
            if strands[i] == strand and min(starts[i] + w, hi) - max(starts[i], lo) >= min_overlap:
                return True
            i += 1
        return False

def source_sites(ccre_path, classes, hits, shift=0):
    """[(site id, ccre id, seq, motif start, strand, score)]: the best strong hit fully inside each cCRE.
    shift moves every cCRE that far first: a control, see cmd_score"""
    sites = []
    with open(ccre_path) as f:
        for line in f:
            p = line.rstrip('\n').split('\t')
            if len(p) < 6 or p[5] not in classes:
                continue
            rec = hits.by_seq.get(p[0])
            if rec is None:
                continue
            starts, strands, scores = rec
            lo, hi = int(p[1]) + shift, int(p[2]) + shift
            i = int(np.searchsorted(starts, lo, 'left'))
            j = int(np.searchsorted(starts, hi - hits.width, 'right'))
            if j > i:
                k = i + int(np.argmax(scores[i:j]))
                sites.append((len(sites), p[4] if len(p) > 4 else p[3], p[0], int(starts[k]), int(strands[k]), int(scores[k])))
    return sites

def score_target(args):
    hal, source, target, sites, hits_t, pad, min_overlap = args
    w = hits_t.width
    with tempfile.NamedTemporaryFile('w', suffix='.bed', delete=False) as bed:
        for sid, ccre, seq, start, strand, score in sites:
            bed.write('{}\t{}\t{}\t{}\n'.format(seq, start, start + w, sid))
        bed_path = bed.name
    psl_path = bed_path + '.psl'
    subprocess.run(['halLiftover', '--outPSLWithName', hal, source, bed_path, target, psl_path],
                   check=True, stderr=subprocess.DEVNULL)
    landings = collections.defaultdict(list)       # site id -> [(aligned bp, seq, start, end, strand)]
    with open(psl_path) as f:
        for line in f:
            p = line.rstrip('\n').split('\t')
            sid, p = int(p[0]), p[1:]
            aligned = sum(int(x) for x in p[18].rstrip(',').split(','))
            landings[sid].append((aligned, p[13], int(p[15]), int(p[16]), -1 if p[8][-1] == '-' else 1))
    os.unlink(bed_path)
    os.unlink(psl_path)

    rows = []
    for sid, ccre, seq, start, strand, score in sites:
        L = landings.get(sid, [])
        if not L:
            rows.append((sid, ccre, target, 0, 0, 0, 0, 0, 0))
            continue
        def conserved(tname, lo, hi, tstrand):
            return hits_t.has_hit(tname, lo - pad, hi + pad, strand * tstrand, min_overlap)
        prim = max(L, key=lambda x: x[0])
        c_prim = conserved(prim[1], prim[2], prim[3], prim[4])
        c_any = c_prim or any(conserved(x[1], x[2], x[3], x[4]) for x in L)
        tlen = hits_t.lengths.get(prim[1], 0)
        bg_hit = bg_n = 0
        for off in CONTROL_OFFSETS:
            lo, hi = prim[2] + off - pad, prim[3] + off + pad
            if lo < 0 or hi > tlen:
                continue
            bg_n += 1
            bg_hit += hits_t.has_hit(prim[1], lo, hi, strand * prim[4], min_overlap)
        rows.append((sid, ccre, target, len(L), prim[0], int(c_prim), int(c_any), bg_hit, bg_n))
    return target, rows

def tally(rows):
    """aligned, conserved at the primary landing, at any landing, background hit and window counts"""
    al = [r for r in rows if r[3] > 0]
    return (len(rows), len(al), sum(r[5] for r in al), sum(r[6] for r in al),
            sum(r[7] for r in al), sum(r[8] for r in al))

def pct(a, b):
    return 100.0 * a / b if b else 0.0

def write_table(out, label_rows):
    out.write('{:28s} {:>8s} {:>16s} {:>17s} {:>9s} {:>8s} {:>10s} {:>9s}\n'.format(
        '', 'sites', 'aligned', 'conserved', 'any-land', 'backgr', 'expected', 'excess'))
    for label, rows in label_rows:
        n, al, cp, ca, bh, bn = tally(rows)
        bg = bh / bn if bn else 0.0
        exp = al * bg
        out.write('{:28s} {:>8,} {:>8,} {:>6.1f}% {:>8,} {:>7.1f}% {:>8.1f}% {:>7.2f}% {:>10,.0f} {:>8.1f}%\n'.format(
            label, n, al, pct(al, n), cp, pct(cp, al), pct(ca, al), 100 * bg, exp, pct(cp - exp, al)))

def cmd_score(args):
    hits_src = Hits(os.path.join(args.hits, args.source + '.npz'), args.pvalue)
    classes = set(args.classes.split(','))
    sites = source_sites(args.ccre, classes, hits_src)
    n_ccre = 0
    with open(args.ccre) as f:
        for line in f:
            p = line.rstrip('\n').split('\t')
            n_ccre += len(p) >= 6 and p[5] in classes
    # The cCREs should have the motif well above chance -- the same windows moved 5 kb do not.  (Only a
    # minority have it at all: V4 calls a cCRE CTCF on CTCF signal in any of hundreds of biosamples.)
    shifted = len(source_sites(args.ccre, classes, hits_src, shift=5000))
    sys.stderr.write('{}: {:,} cCREs of class {}, {:,} with a strong motif ({:.1%}); the same windows moved '
                     '5 kb: {:,} ({:.1%})\n'.format(args.source, n_ccre, args.classes, len(sites),
                                                  len(sites) / max(n_ccre, 1), shifted, shifted / max(n_ccre, 1)))
    if len(sites) < 2 * shifted:
        sys.exit('the cCREs have the motif barely more often than windows 5 kb away: are the cCRE coordinates '
                 'and {} on the same assembly with the same sequence names?'.format(args.source))
    jobs = [(args.hal, args.source, t, sites, Hits(os.path.join(args.hits, t + '.npz'), args.pvalue),
             args.pad, args.min_overlap) for t in args.targets]
    with Pool(min(args.threads, len(jobs))) as pool:
        results = dict(pool.map(score_target, jobs))
    out = sys.stdout
    out.write('# CTCF sites ({} cCREs with a strong MA0139.1 hit, p <= {:g}) from {} lifted through {}\n'.format(
        args.classes, args.pvalue, args.source, args.hal))
    out.write('# aligned: any motif base aligned.  conserved: strong hit at the primary landing on the expected\n'
              '# strand, % of aligned (any-land: at any landing).  backgr: chance of the same in windows 1-10 kb\n'
              '# away.  expected = aligned x backgr; excess = (conserved - expected) as % of aligned.\n')
    write_table(out, [(t, results[t]) for t in args.targets])
    if args.sites:
        with gzip.open(args.sites, 'wt') as f:
            f.write('site\tccre\ttarget\tlandings\taligned_bp\tconserved\tconserved_any\tbg_hits\tbg_windows\n')
            for t in args.targets:
                for r in results[t]:
                    f.write('\t'.join(str(x) for x in r) + '\n')

def cmd_summarize(args):
    data = collections.OrderedDict()
    for path in args.sites:
        # NAME=PATH names the arm NAME; a bare PATH is named by its directory
        arm, eq, path = path.rpartition('=')
        if not eq:
            arm = os.path.basename(os.path.dirname(os.path.abspath(path)))
        per_target = collections.OrderedDict()
        with gzip.open(path, 'rt') as f:
            f.readline()
            for line in f:
                p = line.rstrip('\n').split('\t')
                r = (int(p[0]), p[1], p[2]) + tuple(int(x) for x in p[3:])
                per_target.setdefault(p[2], {})[r[0]] = r
        data[arm] = per_target
    arms = list(data)
    base = args.baseline or arms[0]
    if base not in data:
        sys.exit('baseline {} is not one of {}'.format(base, arms))
    out = sys.stdout
    out.write('# conserved = strong CTCF motif on the expected strand at the landing, % of aligned;\n'
              '# backgr = the chance rate 1-10 kb away; excess = conserved beyond chance, % of aligned\n')
    for t in data[base]:
        out.write('\n## {}\n'.format(t))
        write_table(out, [(a, list(data[a].get(t, {}).values())) for a in arms])
        out.write('-- aligned by the arm but not by {}:\n'.format(base))
        b = data[base].get(t, {})
        write_table(out, [('+ ' + a, [r for k, r in data[a].get(t, {}).items() if r[3] > 0 and b.get(k, (0,) * 4)[3] == 0])
                          for a in arms if a != base])
        out.write('-- aligned by {} but not by the arm:\n'.format(base))
        write_table(out, [('- ' + a, [r for k, r in b.items() if r[3] > 0 and data[a].get(t, {}).get(k, (0,) * 4)[3] == 0])
                          for a in arms if a != base])

def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest='cmd', required=True)
    s = sub.add_parser('scan')
    s.add_argument('--fasta', required=True)
    s.add_argument('--out', required=True)
    s.add_argument('--pvalue', type=float, default=1e-4)
    s = sub.add_parser('score')
    s.add_argument('--hal', required=True)
    s.add_argument('--source', required=True)
    s.add_argument('--ccre', required=True)
    s.add_argument('--hits', required=True)
    s.add_argument('--classes', default='CA-CTCF')
    s.add_argument('--pvalue', type=float, default=1e-4)
    s.add_argument('--pad', type=int, default=10)
    s.add_argument('--min-overlap', type=int, default=10)
    s.add_argument('--threads', type=int, default=4)
    s.add_argument('--sites')
    s.add_argument('targets', nargs='+')
    s = sub.add_parser('summarize')
    s.add_argument('--baseline')
    s.add_argument('sites', nargs='+')
    args = ap.parse_args()
    {'scan': cmd_scan, 'score': cmd_score, 'summarize': cmd_summarize}[args.cmd](args)

if __name__ == '__main__':
    main()
