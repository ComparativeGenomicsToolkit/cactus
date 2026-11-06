"""busco5: progressive cactus on five vertebrates, human to gar, scored three ways against the alignment.

The row's headline numbers come first: cpu_h, busco_recall and busco_precision, and gar_recall and gar_precision,
the same at the root.

  BUSCO     for each ordered pair of genomes (A, B), the BUSCO genes compleasm calls Single in both, lifted base
            by base through the HAL (scoring/busco-align-eval.py).  With the counts of all 20 pairs pooled:
              busco_recall      % of coding bases aligned in phase in the ortholog's CDS
              busco_precision   % of the coding bases aligned at all that are aligned in phase in the ortholog's
                                CDS; the rest land out of phase, outside the CDS, or outside the ortholog's locus
              nc_inlocus        % of aligned intron and flank that lands in the ortholog's locus
            and recall / precision again for the pairs of each split: mouse (human-mouse), opossum (the
            placental-marsupial node), chicken (the amniote node) and gar (the root).  Averages hide most of a
            difference at the deep nodes, so read those.
  CTCF      human CTCF cCREs (ENCODE V4, CA-CTCF) with a strong motif, lifted to each other genome
            (scoring/ctcf-align-eval.py), summed over the four: ctcf_aligned, the sites aligned, and
            ctcf_beyond_chance, the sites whose landing keeps the motif on the right strand beyond the chance
            rate 1-10 kb away.  Alignment of non-homologous sequence raises the first but not the second.
  coverage  halCoverage, a million sampled bases of each genome: the share aligned to each other genome, the mean
            over the 20 pairs.  It rewards aligning noise as much as homology, so it ranks aligners almost exactly
            backwards; it is here to show how much is aligned, not how well
  cost      cpu_h: the alignment's slurm jobs' hours times their cores, summed (retried jobs included);
            align_peak_gib: the most memory any one of them used; cons_h and cons_peak_gib: cactus_consolidated's
            hours summed over the ancestors, and its peak memory; align_h: the alignment's wall-clock hours, mostly
            queueing
  busco_truth  a digest of the BUSCO calls: compleasm does not make the same calls twice, so only rows with the
            same digest are comparable.  They are made once and kept in the cache.

Run to run noise, from repeated arms: BUSCO within about 0.1 point, CTCF within about 15-50 sites, time about 15%.
"""
import collections, glob, gzip, hashlib, os, re, shlex, subprocess

NAME = 'busco5'
DESCRIPTION = 'progressive cactus on human, mouse, opossum, chicken and gar: BUSCO, CTCF and coverage'
# bump when anything below changes the numbers, so rows scored differently are not compared
EVAL_VERSION = 1

SEQFILE = 'busco5.seqfile'
CHROMINFO = 'busco5.chrom-info'
CTCF_SOURCE = 'GCA_000001405.15'   # human: the cCREs are on GRCh38
CCRE_URL = 'https://downloads.wenglab.org/Registry-V4/GRCh38-cCREs.bed'
# --branchScale 2.0 as in the 577-way the tree comes from: its (ROADIES) branch lengths run much shorter than cactus
# expects, so they are doubled.  bar's distances are the scaled ones too.
CACTUS_OPTS = ['--branchScale', '2.0', '--maxCores', '64', '--consCores', '64', '--disableProgress']

COMPLEASM_VERSION = '0.2.9'
LINEAGE = 'vertebrata'
# BUSCO reference proteins used, by NCBI taxon.  compleasm maps every reference protein with miniprot, and
# vertebrata_odb12 has ~1M of them: days per genome.  This panel is ~28k, still covers ~3658 of the 3664
# BUSCOs, and runs hg38 chr21 in under 2 minutes.  vertebrata_odb12 has no human or mouse proteins, so close
# relatives stand in -- primates mouse lemur, marmoset, orangutan; rodents/lagomorph Chinese and golden
# hamster, rabbit -- plus opossum, chicken and spotted gar themselves.
PANEL = ['30608', '9483', '9601', '10029', '10036', '9986', '13616', '9031', '7918']

def check(r):
    for b in ('cactus', 'cactus_consolidated', 'halLiftover', 'halStats', 'halCoverage'):
        if not (os.path.exists(os.path.join(r.bin, b)) or os.path.exists(os.path.join(r.venv, 'bin', b))):
            raise RuntimeError('{} is not in {}: is that install built?'.format(b, r.bin))
    cfg = open(r.path('config.xml')).read()
    # minipoa is a build option; without it bar aborts, but only hours in, once the alignment reaches it
    if re.search(r'\bbaseAligner="minipoa"', cfg) and \
       b'built with minipoa=off' in open(os.path.join(r.bin, 'cactus_consolidated'), 'rb').read():
        raise RuntimeError('baseAligner is minipoa but {}/cactus_consolidated was built with minipoa=off'.format(r.bin))

def run(r):
    from_seqfile = r.seqfile_genomes(r.data(SEQFILE))
    r.step('download', download, r, from_seqfile)
    names = [g for g, _ in from_seqfile]
    r.parallel(lambda: r.step('align', align, r),
               lambda: r.step('compleasm', compleasm_all, r, names),
               lambda: r.step('ctcf-prep', ctcf_prep, r, names))
    r.parallel(lambda: r.step('score-busco', score_busco, r, names),
               lambda: r.step('score-ctcf', score_ctcf, r, names),
               lambda: r.step('score-coverage', score_coverage, r, names))

def local(r):
    """genome -> its downloaded fasta, from the run's seqfile"""
    return collections.OrderedDict(r.seqfile_genomes(r.path(SEQFILE)))

def download(r, from_seqfile):
    """every genome into the cache, and the run's seqfile: the suite's, pointing at the downloads"""
    paths = dict(zip([g for g, _ in from_seqfile], r.fetch_all([u for _, u in from_seqfile])))
    with open(r.data(SEQFILE)) as f, open(r.path(SEQFILE), 'w') as out:
        for line in f:
            p = line.split()
            if len(p) == 2 and p[0] in paths:
                line = '{}\t{}\n'.format(p[0], paths[p[0]])
            out.write(line)

def align(r):
    hal = r.path(NAME + '.hal')
    r.toil('align', ['cactus', r.path(SEQFILE), hal, '--configFile', r.path('config.xml'),
                     '--chromInfo', r.data(CHROMINFO)] + r.cactus_opts(CACTUS_OPTS),
           output=hal)
    if not os.path.exists(hal):
        raise RuntimeError('cactus finished without making {}'.format(hal))

def compleasm_tools(r):
    """compleasm, its own venv (it needs pandas, which a cactus venv need not have) and the lineage, cut down to
    the reference proteins of PANEL: in the cache, once"""
    def setup(d):
        r.sbatch('compleasm-setup', r'''
            cd {d}
            wget -q https://github.com/huangnengCSU/compleasm/releases/download/v{v}/compleasm-{v}_x64-linux.tar.bz2
            tar xjf compleasm-{v}_x64-linux.tar.bz2
            # compleasm unpacks the lineage with one tar.getmember() per HMM, each a scan of the whole gzipped
            # archive: over an hour for vertebrata's 3390 HMMs.  One getmembers() pass does it in seconds.
            python3 - <<'PATCH'
            import pathlib
            p = pathlib.Path('compleasm_kit/compleasm.py'); s = p.read_text()
            old = "tar.extractall(self.download_dir, members=[tar.getmember(u) for u in hmm_files])"
            new = "tar.extractall(self.download_dir, members=[m for m in tar.getmembers() if '.hmm' in m.name])"
            assert s.count(old) == 1, 'compleasm.py has changed; the lineage-unpacking patch no longer applies'
            p.write_text(s.replace(old, new))
            PATCH
            {make_venv}
            venv-compleasm/bin/pip install -q pandas
            venv-compleasm/bin/python compleasm_kit/compleasm.py download {lineage} -L mb_downloads
            db=mb_downloads/{lineage}_odb12/refseq_db.faa
            mv $db.gz $db.full.gz
            gzip -dc $db.full.gz | awk -v panel="{panel}" '
                BEGIN {{ n = split(panel, a, " "); for (i = 1; i <= n; i++) keep[a[i]] = 1 }}
                /^>/ {{ split(substr($1, 2), f, "_"); k = (f[2] in keep) }}
                k' | gzip > $db.gz
            kept=$(gzip -dc $db.gz | grep '>' | sed 's/^>//; s/_.*//' | sort -u | wc -l)
            echo "reference proteins limited to taxa {panel}: $(gzip -dc $db.gz | grep -c '>') covering $kept BUSCOs"
            [ $kept -gt 3500 ] || {{ echo "the taxon panel covers only $kept BUSCOs" >&2; exit 1; }}
            '''.format(d=shlex.quote(d), v=COMPLEASM_VERSION, lineage=LINEAGE, panel=' '.join(PANEL),
                       make_venv=r.make_venv('venv-compleasm')),
            cpus=4, mem='8G', time='4:00:00')
    return r.cached('compleasm-{}-{}-{}'.format(COMPLEASM_VERSION, LINEAGE, '.'.join(PANEL)), setup)

def compleasm_all(r, names):
    """compleasm on each genome, once per download (in the cache), linked into busco/<genome>"""
    tools = compleasm_tools(r)
    fastas = local(r)
    os.makedirs(r.path('busco'), exist_ok=True)
    def one(g):
        def call(d):
            # Each job gets a library of its own -- copies of compleasm's small bookkeeping files, links to the
            # shared lineage -- because every compleasm run re-fetches file_versions.tsv.hash into its library,
            # and runs sharing one read each other's half-written copy, start a re-download and leave a .tmp
            # behind that makes the rest exit ("another process is downloading").
            r.sbatch('compleasm-' + g, '''
                cd {d}
                mkdir lib
                for f in {tools}/mb_downloads/*; do
                    case $f in
                        *.tar.gz) ;;
                        *) if [ -d $f ]; then ln -s $(readlink -f $f) lib/; else cp $f lib/; fi ;;
                    esac
                done
                gzip -dc {fasta} > genome.fa
                {tools}/venv-compleasm/bin/python {tools}/compleasm_kit/compleasm.py run -a genome.fa -o . \\
                    -l {lineage} -L lib -t 32
                rm genome.fa
                '''.format(d=shlex.quote(d), tools=shlex.quote(tools), fasta=shlex.quote(fastas[g]), lineage=LINEAGE),
                cpus=32, mem='64G', time='4:00:00')
        d = r.cached('busco/{}/{}'.format(os.path.basename(tools), digest(fastas[g])), call)
        link(d, r.path('busco', g))
    r.parallel(*[lambda g=g: one(g) for g in names])

def ctcf_prep(r, names):
    """the cCRE registry, and each genome's genome-wide CTCF motif hits (in the cache), linked into ctcf/"""
    os.makedirs(r.path('ctcf', 'hits'), exist_ok=True)
    link(r.fetch(CCRE_URL), r.path('ctcf', 'cCREs.bed'))
    fastas = local(r)
    def one(g):
        def scan(d):
            r.sbatch('ctcf-scan-' + g, 'python3 {} scan --fasta {} --out {}/hits.npz'.format(
                shlex.quote(r.scorer('ctcf-align-eval.py')), shlex.quote(fastas[g]), shlex.quote(d)),
                cpus=1, mem='8G', time='4:00:00')
        d = r.cached('ctcf-scan/v{}/{}'.format(EVAL_VERSION, digest(fastas[g])), scan)
        link(os.path.join(d, 'hits.npz'), r.path('ctcf', 'hits', g + '.npz'))
    r.parallel(*[lambda g=g: one(g) for g in names])

def score_busco(r, names):
    r.sbatch('score-busco', 'python3 {} --hal {} --busco busco --threads 20 {} > busco-eval.txt'.format(
        shlex.quote(r.scorer('busco-align-eval.py')), NAME + '.hal', ' '.join(names)),
        cpus=20, mem='96G', time='4:00:00')

def score_ctcf(r, names):
    targets = [g for g in names if g != CTCF_SOURCE]
    r.sbatch('score-ctcf', 'python3 {} score --hal {} --source {} --ccre ctcf/cCREs.bed --hits ctcf/hits '
             '--sites ctcf-sites.tsv.gz --threads 4 {} > ctcf.txt'.format(
                 shlex.quote(r.scorer('ctcf-align-eval.py')), NAME + '.hal', CTCF_SOURCE, ' '.join(targets)),
             cpus=4, mem='32G', time='4:00:00')

def score_coverage(r, names):
    # The HAL is loaded into memory (--hdf5InMemory): sampling a million random bases out of the HDF5 file on disk
    # is far too slow.  A fixed seed samples the same positions of each genome in every run, so runs are compared
    # on the same bases and the sampling noise cancels out of the comparison.
    os.makedirs(r.path('coverage'), exist_ok=True)
    size = os.path.getsize(r.path(NAME + '.hal'))
    mem = '{}G'.format(int(size * 2 / 2 ** 30) + 8)
    r.parallel(*[lambda g=g: r.sbatch('coverage-' + g, 'halCoverage {} {} --numSamples 1000000 --hdf5InMemory '
                                      '--seed 1 > coverage/{}.tsv'.format(NAME + '.hal', g, g),
                                      cpus=1, mem=mem, time='4:00:00') for g in names])

def digest(path):
    """what a cached result made from the download at path is filed under: the download's own key"""
    return os.path.basename(os.path.dirname(path))

def link(target, name):
    if os.path.lexists(name):
        os.unlink(name)
    os.symlink(target, name)

# --- summary ---

UNALIGNED, ELSEWHERE, LOCUS, OUTPHASE, INPHASE = 0, 1, 2, 3, 4
# a pair of genomes is filed under the deepest split it spans (the most recent common ancestor of the two): human-
# mouse, then each further genome with everything closer to human, so that opossum's pairs show the placental-
# marsupial node, chicken's the amniote node and gar's the root
DEPTH = collections.OrderedDict([('GCA_000001405.15', 0), ('GCA_000001635.9', 0), ('GCF_027887165.1', 1),
                                 ('GCF_016700215.2', 2), ('GCF_040954835.1', 3)])
GROUPS = ['mouse', 'opossum', 'chicken', 'gar']
COMMON = {'GCA_000001405.15': 'human', 'GCA_000001635.9': 'mouse', 'GCF_027887165.1': 'opossum',
          'GCF_016700215.2': 'chicken', 'GCF_040954835.1': 'gar'}

def busco_pairs(path):
    """(A, B) -> {'cds': counts, 'nc': counts} from busco-align-eval.py's #TSV lines, counts indexed as above
    (intron and flank counts summed as nc)"""
    pairs = collections.OrderedDict()
    with open(path) as f:
        for line in f:
            if line.startswith('#TSV\t'):
                _, a, b, n, part, counts, multi = line.rstrip('\n').split('\t')
                d = pairs.setdefault((a, b), {'cds': [0] * 5, 'nc': [0] * 5})
                k = 'cds' if part == 'cds' else 'nc'
                d[k] = [x + int(y) for x, y in zip(d[k], counts.split(','))]
    return pairs

def busco_scores(pairs):
    """(correct, misaligned, nc in-locus) in %, the counts of pairs pooled: correct = coding bases aligned in phase in
    the ortholog's CDS, misaligned = aligned anywhere else (both of all coding bases), nc in-locus = intron and flank
    bases landing in the ortholog's locus, of those aligned"""
    cds, nc = [0] * 5, [0] * 5
    for c in pairs:
        cds = [x + y for x, y in zip(cds, c['cds'])]
        nc = [x + y for x, y in zip(nc, c['nc'])]
    return (pct(cds[INPHASE], sum(cds)), pct(cds[ELSEWHERE] + cds[LOCUS] + cds[OUTPHASE], sum(cds)),
            pct(nc[LOCUS], sum(nc) - nc[UNALIGNED]))

def recall_precision(pairs):
    """(recall, precision) in %, the pairs' counts pooled: coding bases aligned in phase in the ortholog's CDS, of all
    coding bases and of those aligned at all"""
    correct, misaligned, _ = busco_scores(pairs)
    return correct, pct(correct, correct + misaligned)

def busco_groups(pairs):
    """group -> its pairs' counts"""
    out = collections.OrderedDict((g, []) for g in GROUPS)
    for (a, b), c in pairs.items():
        out[GROUPS[max(DEPTH[a], DEPTH[b])]].append(c)
    return out

def busco_truth(r, names):
    """a digest of the BUSCO calls scored against.  compleasm does not give the same calls twice, so scores are only
    comparable between runs with the same calls: the cache's, as long as it is kept"""
    h = hashlib.sha1()
    for g in names:
        tables = [t for t in glob.glob(os.path.join(r.path('busco', g), '**', 'full_table.tsv'), recursive=True)
                  if '/lib/' not in t]
        if len(tables) != 1:
            raise RuntimeError('expected one full_table.tsv under {}, found {}'.format(r.path('busco', g), tables))
        with open(tables[0], 'rb') as f:
            h.update(f.read())
    return h.hexdigest()[:8]

def pct(a, b):
    return 100.0 * a / b if b else 0.0

def ctcf_targets(path):
    """target -> [sites, aligned, conserved at the primary landing, background hits, background windows], tallied
    from ctcf-align-eval.py score's per-site table as its own report tallies them"""
    out = collections.OrderedDict()
    with gzip.open(path, 'rt') as f:
        f.readline()
        for line in f:
            site, ccre, target, landings, aligned_bp, conserved, conserved_any, bg_hits, bg_windows = line.split('\t')
            t = out.setdefault(target, [0, 0, 0, 0, 0])
            t[0] += 1
            if int(landings) > 0:
                t[1] += 1
                t[2] += int(conserved)
                t[3] += int(bg_hits)
                t[4] += int(bg_windows)
    return out

def beyond_chance(t):
    """sites conserved beyond what chance gives at the same number of landings, chance being the motif rate in
    windows 1-10 kb from them"""
    sites, al, cons, bg_hits, bg_windows = t
    return cons - (al * bg_hits / bg_windows if bg_windows else 0.0)

def coverage(r, names):
    """(reference, other) -> % of the reference's sampled bases aligned to the other, from halCoverage"""
    out = collections.OrderedDict()
    for ref in names:
        with open(r.path('coverage', ref + '.tsv')) as f:
            rows = [re.split(r',\s*', l.strip()) for l in f if l.strip()]
        for row in rows[1:]:
            if row[0] != ref and row[0] in names:
                out[(ref, row[0])] = pct(int(row[1]), 1000000)
    return out

def consolidated(r):
    """ancestor -> (hours, peak GiB) of its last successful cactus_consolidated, from the cactus log: the cost of the
    alignment itself, which the run's wall-clock time, mostly queueing, does not show"""
    per = collections.OrderedDict()
    for log in (r.path('logs', 'align.out'), r.path('logs', 'align.log')):
        if os.path.exists(log):
            with open(log, errors='replace') as f:
                for line in f:
                    m = re.search(r'Successfully ran cactus_consolidated\((\w+)\).* in ([0-9.]+) seconds and ([0-9.]+) ?([KMGT])i',
                                  line)
                    if m:
                        per[m.group(1)] = (float(m.group(2)) / 3600,
                                           float(m.group(3)) * {'K': 2 ** -20, 'M': 2 ** -10, 'G': 1, 'T': 2 ** 10}[m.group(4)])
            if per:
                break
    return per

def summarize(r):
    names = list(DEPTH)
    pairs = busco_pairs(r.path('busco-eval.txt'))
    if sorted(pairs) != sorted((a, b) for a in names for b in names if a != b):
        raise RuntimeError('busco-eval.txt does not have every pair of {}'.format(', '.join(names)))
    groups = busco_groups(pairs)
    root = GROUPS[-1]
    m = collections.OrderedDict()
    cpu_h, peak_gib = r.slurm_cost('align')
    m['cpu_h'] = cpu_h
    m['busco_recall'], m['busco_precision'] = recall_precision(pairs.values())
    m[root + '_recall'], m[root + '_precision'] = recall_precision(groups[root])
    for g in GROUPS[:-1]:
        m[g + '_recall'], m[g + '_precision'] = recall_precision(groups[g])
    m['nc_inlocus'] = busco_scores(pairs.values())[2]
    ctcf = ctcf_targets(r.path('ctcf-sites.tsv.gz'))
    m['ctcf_aligned'] = sum(t[1] for t in ctcf.values())
    m['ctcf_beyond_chance'] = int(round(sum(beyond_chance(t) for t in ctcf.values())))
    cov = coverage(r, names)
    m['coverage'] = sum(cov.values()) / len(cov)
    for k in m:
        if isinstance(m[k], float):
            m[k] = round(m[k], 2)
    cons = consolidated(r)
    m['cons_h'] = round(sum(h for h, gib in cons.values()), 2) if cons else ''
    m['cons_peak_gib'] = round(max(gib for h, gib in cons.values()), 1) if cons else ''
    m['align_h'] = r.hours('align')
    m['align_peak_gib'] = peak_gib
    m['busco_truth'] = busco_truth(r, names)
    return m

def compare(runs, out):
    """the BUSCO groups and cactus_consolidated per ancestor, a column per run; then each scorer's own side-by-side
    tables: BUSCO per pair, coverage per pair, and CTCF per target with what each run aligns that the first does not"""
    names = list(DEPTH)
    labels = runs[0].short_labels([r.p['label'] for r in runs])
    head = lambda title: out.write('{:24s}'.format(title) + ''.join(' {:>17s}'.format(l[:17]) for l in labels) + '\n')
    out.write('# BUSCO, pairs pooled by the deepest split they span: % of CDS correct / misaligned\n')
    head('group')
    groups = [busco_groups(busco_pairs(r.path('busco-eval.txt'))) for r in runs]
    for g in GROUPS:
        out.write('{:24s}'.format(g) + ''.join(' {:>8.2f} / {:>5.2f}'.format(*busco_scores(gr[g])[:2]) for gr in groups) + '\n')
    cons = [consolidated(r) for r in runs]
    out.write('\n# cactus_consolidated per ancestor: hours / peak GiB\n')
    head('ancestor')
    for anc in sorted(set().union(*cons)):
        out.write('{:24s}'.format(anc) + ''.join(' {:>8.2f} / {:>5.0f}'.format(*c[anc]) if anc in c else ' {:>17s}'.format('-')
                                                 for c in cons) + '\n')
    covs = [coverage(r, names) for r in runs]
    out.write('\n# coverage, per pair: % of the reference\'s bases aligned to the other\n')
    head('reference -> other')
    for pair in covs[0]:
        out.write('{:24s}'.format('{} -> {}'.format(*[COMMON[g] for g in pair])) +
                  ''.join(' {:>16.2f}%'.format(c[pair]) if pair in c else ' {:>17s}'.format('-') for c in covs) + '\n')
    out.write('\n')
    out.flush()
    py = os.path.join(runs[0].venv, 'bin', 'python')
    subprocess.run([py, runs[0].scorer('busco-align-eval.py'), '--summarize'] +
                   ['{}={}'.format(l, r.path('busco-eval.txt')) for l, r in zip(labels, runs)], stdout=out, check=True)
    out.write('\n')
    out.flush()
    subprocess.run([py, runs[0].scorer('ctcf-align-eval.py'), 'summarize', '--baseline', labels[0]] +
                   ['{}={}'.format(l, r.path('ctcf-sites.tsv.gz')) for l, r in zip(labels, runs)], stdout=out, check=True)
