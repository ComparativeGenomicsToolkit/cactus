"""hprc30: cactus-pangenome on 30 haplotypes -- CHM13, GRCh38 and both haplotypes of 14 HPRC release-2 samples, HG002
among them, as in the HPRC v2.1 graph -- with CHM13 as reference, scored by the accuracy of its VCFs.

  accuracy  toil-vcfcomp (aardvark, genotype comparison, no normalization) on two of the graph's VCFs:
              raw10k  the raw VCF less its sites over 10 kb (vcfbub -l 0 -a 10000), since aardvark skips any allele
                      that long, so the default VCF's -r 100000 would bias it.  The headline: vcfwave hides under-
                      alignment in the graph
              wave    vcfwave's VCF
            against two truths, both in-sample (the samples are in the graph):
              dipcall  each sample's own dipcall calls on CHM13v2 (HPRC), within dipcall's confident regions.
                       dipcall_f1 is the mean GT F1 over the 14 samples; dipcall_err the errors (truth FN + query FP)
                       summed over them, which move more than the mean F1 does; dipcall_sv50_err the same for
                       variants of 50 bp and over
              GIAB     HG002 against GIAB's Q100 benchmark on CHM13v2.0: small variants (smvar) and structural
                       variants (stvar), GT F1 and errors
            Run to run noise is about 0.0002 F1.
  graph     nodes, edges and Mb of the graph and the % of input sequence left out of it (cactus-graphmap-join's
            stats; blank for a cactus that does not write them)
  cost      build_h: cactus-pangenome's wall-clock hours, queueing included

The truths are fixed manifests (hprc30.dipcall.csv, hprc30.giab.csv) and toil-vcfcomp a fixed commit, so every run
is scored the same way and runs compare by their rows.  A run takes about 130 GB, almost all of it in hprc30/:
delete its chrom-* directories once it is recorded, unless they are wanted.  cactus builds that predate per-job
walltimes (--doubleTime) need --toilOpts '--slurmTime 24:00:00' or so: Toil's default is too short for some jobs.
"""
import collections, csv, os, re, shlex
import xml.etree.ElementTree as ET

NAME = 'hprc30'
DESCRIPTION = 'cactus-pangenome on CHM13, GRCh38 and 28 HPRC haplotypes: dipcall and GIAB HG002 accuracy'
# bump when anything below changes the numbers, so rows scored differently are not compared
EVAL_VERSION = 1

SEQFILE = 'hprc30.seqfile'
# truthset name in the manifest -> manifest
TRUTH = collections.OrderedDict([('dipcall', 'hprc30.dipcall.csv'), ('giab', 'hprc30.giab.csv')])
REF_FASTA = 'https://hgdownload.soe.ucsc.edu/goldenPath/hs1/bigZips/hs1.fa.gz'   # CHM13v2, as the truths
TOIL_VCFCOMP_REPO = 'https://github.com/glennhickey/toil-vcfcomp.git'
TOIL_VCFCOMP_COMMIT = '977f349326946b6a2a8483b350a8096c790b9233'
CACTUS_OPTS = ['--reference', 'CHM13', 'GRCh38', '--mgSplit', '--gbz', '--vcf', '--vcfwave', '--snarlStats',
               '--disableProgress', '--binariesMode', 'local', '--consCores', '32', '--mapCores', '12',
               '--mgCores', '12', '--indexCores', '32', '--vcfwaveCores', '18']
QUERIES = ('raw10k', 'wave')
# The maternal haplotypes the HPRC v2.1 graph was built from are v1.2.0: v1.0.1 (public, in the seqfile) with the
# mitochondrion replaced by GenBank's corrected record, ACCESSION.2 for ACCESSION.1.  v1.2.0 itself is not public, so
# it is remade: genome -> (accession, md5 of the .2 record's sequence, uppercased), the md5 checked against the
# mitochondria of a 30-way built from the v1.2.0 files.
MITO = collections.OrderedDict([
    ('HG00408.2', ('JBHDVL010000072', '90c408e792379f08fc3d137cd7fb6089')),
    ('HG00597.2', ('JBHDUV010000091', 'cdeff7e6df333d88d411fe5d339dbd1e')),
    ('HG01192.2', ('JBHDVJ010000069', 'aaf7e099ea4affde5b0eeceb169210a6')),
    ('HG01261.2', ('JBHDVH010000084', 'c0d92ac63023050a5b1f81d4efb7ca28')),
    ('HG02015.2', ('JBHDVF010000101', '08821a074c0b7ed66732a44afdf88498')),
    ('HG02056.2', ('JBHDVD010000077', 'ecb324a2417dc4b923ce0a8c60eed9a0')),
    ('HG02129.2', ('JBHDVB010000071', '4dee8e545c8452c26a5d6b0e83ca80a1')),
    ('HG01975.2', ('JBJUVT010000092', 'b57786e28f2384a06f52181efdf00a87')),
    ('HG02602.2', ('JBJUVL010000063', '5fdb5632e2ae35b8db5193a012449022')),
    ('HG04187.2', ('JBJUVP010000104', '5f54d1a6b9d0ad78a2188a67b10a1891')),
    ('NA18879.2', ('JBJUVN010000087', 'e317fd44ed42065f0a4c6b939e0dbd78')),
    ('NA20752.2', ('JBJUVR010000089', 'e75509379b6fd1f2a21e01bd974709bb')),
    ('HG00232.2', ('JBHIKN010000065', 'bc09d683608f5039b8a758bfa834b97d'))])
EFETCH = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id={}&rettype=fasta&retmode=text'

def prefix(r):
    return r.path(NAME, NAME)

def check(r):
    for b in ('cactus-pangenome', 'minigraph', 'vg', 'vcfwave', 'vcfbub', 'bgzip', 'tabix'):
        if not (os.path.exists(os.path.join(r.bin, b)) or os.path.exists(os.path.join(r.venv, 'bin', b))):
            raise RuntimeError('{} is not in {}: is that install built?'.format(b, r.bin))

def run(r):
    r.step('download', download, r)
    r.parallel(lambda: r.step('build', build, r),
               lambda: r.step('toil-vcfcomp', toil_vcfcomp, r))
    r.step('raw10k', raw10k, r)
    r.parallel(*[lambda t=t: r.step('score-' + t, score, r, t) for t in TRUTH])

def download(r):
    """every assembly into the cache, the v1.2.0 haplotypes remade, and the run's seqfile: the suite's, pointing at
    them"""
    genomes = r.seqfile_genomes(r.data(SEQFILE))
    paths = dict(zip([g for g, _ in genomes], r.fetch_all([u for _, u in genomes])))
    r.parallel(*[lambda g=g: paths.__setitem__(g, mito_corrected(r, g, paths[g])) for g in MITO])
    pansn = reads_pansn(r)
    with open(r.data(SEQFILE)) as f, open(r.path(SEQFILE), 'w') as out:
        for line in f:
            p = line.split()
            if len(p) == 2 and p[0] in paths:
                line = '{}\t{}\n'.format(re.sub(r'\.(\d+)$', r'#\1', p[0]) if pansn else p[0], paths[p[0]])
            out.write(line)

def reads_pansn(r):
    """whether the cactus under test takes a haplotype as SAMPLE#N, the PanSN spelling, in a seqfile.  The suite names
    them SAMPLE.N, the only spelling a cactus older than <graphmap dotHaplotypes> reads as a haplotype; one with it
    reads SAMPLE#N, and is to stop reading SAMPLE.N (dotHaplotypes="0" is to be the default)"""
    return any('dotHaplotypes' in e.attrib for e in ET.parse(r.path('config.xml')).iter('graphmap'))

def mito_corrected(r, genome, v101):
    """genome's v1.2.0 fasta, made from its v1.0.1 download (in the cache, once)"""
    sample, hap = genome.split('.')
    acc, md5 = MITO[genome]
    name = os.path.basename(v101).replace('_v1.0.1.', '_v1.2.0.')
    def make(d):
        r.sbatch('mito-' + genome, '''
            cd {d}
            # the 13 jobs start together, and NCBI turns away more than 3 requests a second
            wget -nv --tries 10 --waitretry 60 --retry-on-http-error=429,500,502,503,504 -O mito.fa '{url}'
            [ "$(head -1 mito.fa | cut -c2- | cut -d' ' -f1)" = {acc}.2 ] || {{ echo "GenBank did not return {acc}.2" >&2; exit 1; }}
            seq=$(grep -v '^>' mito.fa | tr -d '\\n' | tr a-z A-Z)
            [ "$(printf %s "$seq" | md5sum | cut -d' ' -f1)" = {md5} ] || {{ echo "{acc}.2 is not the sequence HPRC v2.1 has" >&2; exit 1; }}
            gzip -dc {src} | awk -v drop='{old}' '/^>/ {{ skip = (substr($1, 2) == drop); n += skip }} !skip
                END {{ if (n != 1) {{ print "found " n + 0 " contigs named " drop > "/dev/stderr"; exit 1 }} }}' > genome.fa
            printf '>%s\\n%s\\n' '{new}' "$seq" >> genome.fa
            bgzip -@ 4 genome.fa
            mv genome.fa.gz {name}
            rm mito.fa
            '''.format(d=shlex.quote(d), url=EFETCH.format(acc + '.2'), acc=acc, md5=md5, src=shlex.quote(v101),
                       old='{}#{}#{}.1'.format(sample, hap, acc), new='{}#{}#{}.2'.format(sample, hap, acc),
                       name=shlex.quote(name)), cpus=4, mem='4G', time='2:00:00')
    return os.path.join(r.cached('hprc-v1.2.0/{}-{}'.format(genome, os.path.basename(os.path.dirname(v101))), make), name)

def build(r):
    os.makedirs(r.path(NAME), exist_ok=True)
    r.toil('build', ['cactus-pangenome', r.path(SEQFILE), '--outDir', r.path(NAME), '--outName', NAME,
                     '--configFile', r.path('config.xml')] + CACTUS_OPTS +
           shlex.split(r.p.get('cactus_opts') or ''), output=prefix(r) + '.gbz')
    for q in ('raw', 'wave'):
        if not os.path.exists('{}.{}.vcf.gz'.format(prefix(r), q)):
            raise RuntimeError('cactus-pangenome finished without making {}.{}.vcf.gz'.format(prefix(r), q))

def raw10k(r):
    # -a, not -r: -r misses ALT-only alleles over 10 kb.  vcfbub cannot read stdin.
    r.sbatch('raw10k', '''
        vcfbub -l 0 -a 10000 -i {p}.raw.vcf.gz | bgzip -@ 4 > {p}.raw10k.vcf.gz.tmp
        mv {p}.raw10k.vcf.gz.tmp {p}.raw10k.vcf.gz
        tabix -p vcf {p}.raw10k.vcf.gz
        '''.format(p=shlex.quote(prefix(r))), cpus=4, mem='8G', time='2:00:00')

def toil_vcfcomp(r):
    """toil-vcfcomp at TOIL_VCFCOMP_COMMIT in a venv of its own (on the shared filesystem: Toil does not ship it to
    the workers), in the cache, once"""
    def install(d):
        r.sh(['bash', '-c', '''set -eo pipefail
            git clone -q {repo} {d}/src
            git -C {d}/src checkout -q {commit}
            {make_venv}
            . {d}/venv/bin/activate
            pip install -q -U pip
            pip install -q -U -r {d}/src/toil-requirement.txt
            pip install -q -U {d}/src
            which toil-vcfcomp _toil_worker
            '''.format(repo=TOIL_VCFCOMP_REPO, d=shlex.quote(d), commit=TOIL_VCFCOMP_COMMIT,
                       make_venv=r.make_venv(os.path.join(d, 'venv')))], 'toil-vcfcomp')
    return os.path.join(r.cached('toil-vcfcomp-' + TOIL_VCFCOMP_COMMIT[:7], install), 'venv')

def score(r, truth):
    venv = toil_vcfcomp(r)
    # The workers need a bcftools with the +split plugin, and cactus's has none: every cactus install comes off the
    # PATH, and a compute node's bcftools is checked before anything is launched.
    r.sbatch('bcftools-' + truth, 'type -a bcftools; bcftools plugin -l 2>&1 | grep -w split',
             cpus=1, mem='1G', time='0:10:00', venv=venv, cactus_bin=False)
    out = r.path('vcfcomp-' + truth)
    r.toil('score-' + truth, ['toil-vcfcomp', '--queryVcf'] +
           ['{}={}.{}.vcf.gz'.format(q, prefix(r), q) for q in QUERIES] +
           ['--truthIndex', r.data(TRUTH[truth]), '--ref', REF_FASTA, '--outDir', out, '--evalType', 'aardvark',
            '--noNormalize', '--binariesMode', 'singularity', '--threads', '8', '--compareMemory', '32Gi',
            '--maxJobs', '300', '--slurmTime', '2:00:00'],
           output=os.path.join(out, 'report.tsv'), env=r.env(venv, cactus_bin=False))
    with open(os.path.join(out, 'report.txt')) as f:
        if not any(l.startswith('normalize : no') for l in f):
            raise RuntimeError('{}/report.txt does not say "normalize : no"'.format(out))

# --- summary ---

def report(r, truth):
    with open(r.path('vcfcomp-' + truth, 'report.tsv')) as f:
        return list(csv.DictReader(f, delimiter='\t'))

def gt(rows, callset, truthset, variant_type='ALL'):
    """per sample, aardvark's GT comparison of callset against truthset, over all regions"""
    return collections.OrderedDict((x['sample'], x) for x in rows if
        x['callset'] == callset and x['truthset'] == truthset and x['tool'] == 'aardvark' and
        x['comparison'] == 'GT' and x['region_label'] == 'ALL' and x['filter'] == 'ALL' and
        x['variant_type'] == variant_type)

def errors(x):
    return int(x['truth_fn']) + int(x['query_fp'])

def f1_err(per_sample):
    """(mean F1, summed errors)"""
    xs = list(per_sample.values())
    return sum(float(x['metric_f1']) for x in xs) / len(xs), sum(errors(x) for x in xs)

def chrom(r, truth, callset, metrics):
    """contig -> count of metrics summed over the samples, from report.chrom.tsv"""
    out = collections.Counter()
    with open(r.path('vcfcomp-' + truth, 'report.chrom.tsv')) as f:
        for x in csv.DictReader(f, delimiter='\t'):
            if x['callset'] == callset and x['truthset'] == truth and x['tool'] == 'aardvark' and x['metric'] in metrics:
                out[x['contig']] += int(x['count'])
    return out

def graph(r):
    """(nodes, edges, Mb, % of input left out) from cactus-graphmap-join's stats, if it wrote them"""
    stats = prefix(r) + '.stats'
    n = e = l = 0
    try:
        with open(os.path.join(stats, 'graph-stats.tsv')) as f:
            for row in csv.reader(f, delimiter='\t'):
                if row and not row[0].startswith('#'):
                    n, e, l = n + int(row[1]), e + int(row[2]), l + int(row[3])
        with open(os.path.join(stats, 'clipped-by-genome.tsv')) as f:
            clipped = [float(row[4]) for row in csv.reader(f, delimiter='\t') if row[:2] == ['TOTAL', 'TOTAL']]
    except FileNotFoundError:
        return '', '', '', ''
    return n, e, round(l / 1e6, 1), clipped[0] if clipped else ''

def summarize(r):
    m = collections.OrderedDict()
    dipcall, giab = report(r, 'dipcall'), report(r, 'giab')
    samples = len(open(r.data(TRUTH['dipcall'])).readlines()) - 1
    for q in QUERIES:
        per = gt(dipcall, q, 'dipcall')
        if len(per) != samples:
            raise RuntimeError('{} has {} samples scored against dipcall, not {}'.format(q, len(per), samples))
        m['dipcall_f1_' + q], m['dipcall_err_' + q] = f1_err(per)
        if q == 'raw10k':
            m['dipcall_sv50_err_' + q] = chrom(r, 'dipcall', q, ('SV50-FN', 'SV50-FP'))['_total_']
    for q in QUERIES:
        for t, label in (('v5.0q_smvar', 'smvar'), ('v5.0q_stvar', 'stvar')):
            per = gt(giab, q, t)
            if list(per) != ['HG002']:
                raise RuntimeError('{} against {}: expected HG002 alone, got {}'.format(q, t, list(per)))
            m['{}_f1_{}'.format(label, q)], err = f1_err(per)
            if q == 'raw10k':
                m['{}_err_{}'.format(label, q)] = err
    for k in m:
        if isinstance(m[k], float):
            m[k] = round(m[k], 6)
    m['nodes'], m['edges'], m['graph_mb'], m['clipped_pct'] = graph(r)
    m['build_h'] = r.hours('build')
    return m

def contig_order(c):
    """_total_, then chr1-22, then the rest by name"""
    m = re.match(r'(?:.*#)?chr(\d+)$', c)
    return (c != '_total_', 0 if m else 1, int(m.group(1)) if m else 0, c)

def compare(runs, out):
    """dipcall errors per sample and per chromosome, GIAB F1 per variant type: a column per run, then the change
    from the first run"""
    labels = [r.p['label'] for r in runs]
    def head(title):
        out.write('{:24s}'.format(title) + ''.join(' {:>14s}'.format(l[:14]) for l in labels) +
                  ''.join(' {:>9s}'.format('\u0394' + l[:8]) for l in labels[1:]) + '\n')
    def row(name, vals, fmt, pct=True):
        base = vals[0]
        out.write('{:24s}'.format(name) + ''.join(' {:>14}'.format(fmt.format(v) if v is not None else '-') for v in vals) +
                  ''.join(' {:>9s}'.format('-' if v is None or base is None else
                          '{:+.2f}%'.format(100.0 * (v - base) / base) if pct and base else '{:+.6f}'.format(v - base))
                          for v in vals[1:]) + '\n')
    dip = [report(r, 'dipcall') for r in runs]
    for q in QUERIES:
        per = [gt(d, q, 'dipcall') for d in dip]
        out.write('# dipcall, {}: errors (truth FN + query FP) per sample\n'.format(q))
        head('sample')
        for s in per[0]:
            row(s, [errors(p[s]) if s in p else None for p in per], '{:,}')
        row('total', [sum(errors(x) for x in p.values()) for p in per], '{:,}')
        for i, p in enumerate(per[1:], 1):
            better = sum(1 for s in per[0] if s in p and errors(p[s]) < errors(per[0][s]))
            worse = sum(1 for s in per[0] if s in p and errors(p[s]) > errors(per[0][s]))
            out.write('{}: fewer errors than {} in {} of {} samples, more in {}\n'.format(
                labels[i], labels[0], better, len(per[0]), worse))
        out.write('\n')
    out.write('# dipcall, raw10k: errors per chromosome, summed over the samples\n')
    head('contig')
    per = [chrom(r, 'dipcall', 'raw10k', ('ERRORS',)) for r in runs]
    for c in sorted(per[0], key=contig_order):
        row(c, [p.get(c) for p in per], '{:,}')
    giab = [report(r, 'giab') for r in runs]
    for t in ('v5.0q_smvar', 'v5.0q_stvar'):
        out.write('\n# GIAB HG002 {}, raw10k: GT F1 by variant type\n'.format(t))
        head('variant_type')
        types = [x['variant_type'] for x in giab[0] if x['callset'] == 'raw10k' and x['truthset'] == t and
                 x['tool'] == 'aardvark' and x['comparison'] == 'GT' and x['region_label'] == 'ALL' and x['filter'] == 'ALL']
        for vt in types:
            vals = [gt(g, 'raw10k', t, vt).get('HG002') for g in giab]
            if not any(v and v['metric_f1'] for v in vals):
                continue
            row(vt, [float(v['metric_f1']) if v and v['metric_f1'] else None for v in vals], '{:.6f}', pct=False)
