#!/usr/bin/env python3
"""Big evaluations of cactus on a slurm cluster, from public data only, summarized into one table row per run.

A suite (suites/<name>.py) aligns a fixed set of public genomes with the cactus install under test, scores
the alignment, and boils the scores down to a row of numbers.  One run is one install, with one config, on
one suite; compare installs or configs by running each and comparing their rows (record, compare).  Runs
share a cache of everything that does not depend on the install: downloads, BUSCO calls, motif scans and
the scoring tools.

  bigeval.py list
  bigeval.py run SUITE --cactus VENV --out DIR [--label NAME] [--set [ELEMENT.]ATTR=VALUE ...] [--detach]
  bigeval.py resume DIR [--detach]
  bigeval.py status DIR
  bigeval.py stop DIR
  bigeval.py summarize DIR
  bigeval.py record DIR [--results TSV]
  bigeval.py compare DIR DIR ...

VENV is a cactus virtualenv inside its install, the layout both a source build and a binary release have:
its binaries are in the bin/ beside it (cactus.beta/venv-cactus-beta -> cactus.beta/bin).  Run this on a
head node that can submit to slurm and reach the internet: the driver and Toil leaders run here, and
everything else goes to slurm, so nothing here needs AVX2 or much memory.  A run takes days; --detach
puts it in the background (its own session, logging to DIR/driver.log).  A run that stopped picks up
where it left off with resume: finished steps are skipped and Toil workflows restarted.

Results: DIR/summary.tsv (the row) and DIR/summary.json (the row and everything behind it).  record
appends the row to results/<suite>.tsv beside this script, which is tracked in git.
"""
import argparse, collections, contextlib, csv, datetime, fcntl, hashlib, importlib.util, json, os, re, shlex
import shutil, signal, socket, subprocess, sys, textwrap, threading, time, traceback, urllib.parse, urllib.request
import xml.etree.ElementTree as ET
from concurrent.futures import ThreadPoolExecutor

HARNESS = os.path.dirname(os.path.abspath(__file__))
SUITES = os.path.join(HARNESS, 'suites')
SCORING = os.path.join(HARNESS, 'scoring')
RESULTS = os.path.join(HARNESS, 'results')
# columns every row starts with, before the suite's metrics
ID_COLUMNS = ['date', 'label', 'cactus_commit', 'cactus_version', 'dirty', 'config', 'cactus_opts', 'eval']

def now():
    return datetime.datetime.now().astimezone().isoformat(timespec='seconds')

def hms(seconds):
    seconds = int(seconds)
    return '{}:{:02d}:{:02d}'.format(seconds // 3600, seconds % 3600 // 60, seconds % 60)

def read_json(path):
    try:
        with open(path) as f:
            return json.load(f, object_pairs_hook=collections.OrderedDict)
    except FileNotFoundError:
        return None

def write_json(path, obj):
    with open(path + '.tmp', 'w') as f:
        json.dump(obj, f, indent=2)
        f.write('\n')
    os.replace(path + '.tmp', path)

def suite_names():
    return sorted(f[:-3] for f in os.listdir(SUITES) if f.endswith('.py') and not f.startswith('_'))

def load_suite(name):
    path = os.path.join(SUITES, name + '.py')
    if not os.path.exists(path):
        sys.exit('no suite {}; there are: {}'.format(name, ', '.join(suite_names())))
    spec = importlib.util.spec_from_file_location('bigeval_suite_' + name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod

def venv_python(venv, code):
    return subprocess.run([os.path.join(venv, 'bin', 'python'), '-c', code], check=True,
                          capture_output=True, text=True).stdout.strip()

def git(d, *args):
    r = subprocess.run(['git', '-C', d] + list(args), capture_output=True, text=True)
    return r.stdout.strip() if r.returncode == 0 else None

def identity(venv):
    """what is being tested: the installed cactus package's commit, and the checkout it sits in if it is one.
    A venv installed without -e keeps the code it was installed with whatever the checkout does since, so the
    installed commit is the one that counts; a checkout at another commit, or with local changes, is noted."""
    ident = collections.OrderedDict()
    try:
        ident['cactus_commit'] = venv_python(venv, 'from cactus.shared.version import cactus_commit; print(cactus_commit)')
    except subprocess.CalledProcessError:
        ident['cactus_commit'] = 'unknown'
    checkout = os.path.dirname(os.path.abspath(venv))
    ident['checkout'] = checkout
    head = git(checkout, 'rev-parse', 'HEAD')
    ident['checkout_head'] = head or ''
    ident['cactus_version'] = (git(checkout, 'describe', '--tags', '--always') or '') if head else ''
    # changes to the code, not those its build makes (cactus's Makefile edits some submodules' Makefiles; a submodule
    # at another commit still counts) nor a binary release's missing files (it ships part of the tree, with its .git)
    ident['dirty'] = bool(head and git(checkout, 'diff', 'HEAD', '--name-only', '--diff-filter=d',
                                       '--ignore-submodules=dirty'))
    return ident

def harness_identity():
    head = git(HARNESS, 'rev-parse', 'HEAD')
    return {'harness_commit': head or 'unknown',
            'harness_dirty': bool(head and git(HARNESS, 'status', '--porcelain', '-uno', '--', '.'))}

def stock_config(venv):
    return venv_python(venv, 'import os, cactus; print(os.path.join(os.path.dirname(cactus.__file__), '
                             '"cactus_progressive_config.xml"))')

def edit_config(src, sets, dst):
    """copy cactus config src to dst with each [ELEMENT.]ATTR=VALUE of sets applied.  The attribute must already
    be in the config, on exactly one element (of tag ELEMENT, if given): a typo, or a build that predates the
    option, fails here rather than as a run that quietly ignored it"""
    with open(src) as f:
        text = f.read()
    for s in sets:
        key, eq, value = s.partition('=')
        tag, _, attr = key.rpartition('.')
        if not eq or not attr:
            raise ValueError('--set {}: expected [ELEMENT.]ATTR=VALUE'.format(s))
        hits = [e for e in ET.fromstring(text).iter(tag or None) if attr in e.attrib]
        if len(hits) != 1:
            raise ValueError('--set {}: {} {} elements with a {} attribute in {}{}'.format(
                s, len(hits), '<{}>'.format(tag) if tag else 'config', attr, src,
                '; name the element, e.g. {}.{}'.format(hits[0].tag, attr) if len(hits) > 1 and not tag else ''))
        # Edit the text rather than write the tree back, so the config keeps its layout and comments and diffs
        # cleanly against the one it came from.  Comments are blanked out (keeping offsets) before the start tags
        # are searched, so that an attribute quoted in one is not taken for the real thing.
        masked = re.sub(r'<!--.*?-->', lambda m: ' ' * len(m.group()), text, flags=re.S)
        spans = []
        for m in re.finditer(r'<([A-Za-z_][\w.:-]*)\s[^<>]*>', masked):
            a = re.search(r'\s{}\s*=\s*(["\'])(.*?)\1'.format(re.escape(attr)), m.group(0), re.S)
            if a and (not tag or m.group(1) == tag):
                spans.append((m.start() + a.start(2), m.start() + a.end(2), a.group(1)))
        if len(spans) != 1:
            raise ValueError('--set {}: found the attribute but could not edit it in {}'.format(s, src))
        start, end, quote = spans[0]
        escaped = value.replace('&', '&amp;').replace('<', '&lt;').replace(quote, '&quot;' if quote == '"' else '&apos;')
        text = text[:start] + escaped + text[end:]
        if [e.get(attr) for e in ET.fromstring(text).iter(tag or None) if attr in e.attrib] != [value]:
            raise ValueError('--set {}: the edit did not take in {}'.format(s, src))
    with open(dst, 'w') as f:
        f.write(text)

def seqfile_genomes(path):
    """(name, path or URL) of each genome of a cactus seqfile: tree first, then a line per genome"""
    genomes = []
    with open(path) as f:
        for line in f:
            p = line.split()
            if len(p) == 2 and not p[0].startswith('(') and not p[0].startswith('#'):
                genomes.append((p[0].lstrip('*'), p[1]))
    return genomes

def make_venv(d):
    """shell that makes a virtualenv with pip at d, with the python3 on the PATH: in a step, the cactus venv's.  A
    python without ensurepip (Debian's, short of python3-venv) cannot put pip in a venv itself, so the cactus venv's
    pip does it"""
    return ('python3 -m venv {d} || {{ rm -rf {d}; python3 -m venv --without-pip {d} && '
            'python3 -m pip -q --python {d}/bin/python install pip; }}').format(d=shlex.quote(d))

def slurm_seconds(t):
    """a slurm time limit in seconds: MM, MM:SS, HH:MM:SS or D-HH[:MM[:SS]], or infinite"""
    if t in ('infinite', 'UNLIMITED'):
        return float('inf')
    days, dash, rest = t.rpartition('-')
    f = [int(x) for x in rest.split(':')]
    if dash:
        h, m, s = f + [0] * (3 - len(f))
    else:
        h, m, s = (0, f[0], 0) if len(f) == 1 else [0] * (3 - len(f)) + f
    return ((int(days or 0) * 24 + h) * 60 + m) * 60 + s

class Run:
    """one run's directory, parameters and the means to run its steps:

      DIR/run.json        the parameters, and what was tested (identity)
      DIR/config.xml      the cactus config it aligns with
      DIR/steps/S.json    each step's status and times; a done step is not run again
      DIR/logs/           a log per step (Toil's --logFile as S.log, everything else as S.out)
      DIR/js/             Toil jobstores
      DIR/jobs/           the scripts given to sbatch
    """
    def __init__(self, out, params=None):
        self.out = os.path.abspath(out)
        self.p = params or read_json(self.path('run.json'))
        if self.p is None:
            sys.exit('{} is not a bigeval run (no run.json)'.format(self.out))
        self.suite = load_suite(self.p['suite'])
        self.venv = self.p['cactus']
        self.bin = os.path.join(os.path.dirname(self.venv), 'bin')
        self.cache = self.p['cache']
        self.tag = 'bigeval.' + os.path.basename(self.out)
        self._print_lock = threading.Lock()
        self._env_lock = threading.Lock()
        self._activated = {}
        self._slurm_lock = threading.Lock()
        self._partitions = None

    seqfile_genomes = staticmethod(seqfile_genomes)
    make_venv = staticmethod(make_venv)

    def path(self, *p):
        return os.path.join(self.out, *p)

    def cache_path(self, *p):
        return os.path.join(self.cache, *p)

    def data(self, name):
        """a file shipped with the suites"""
        return os.path.join(SUITES, name)

    def scorer(self, name):
        return os.path.join(SCORING, name)

    def log(self, msg):
        with self._print_lock:
            print('[{}] {}'.format(now(), msg), flush=True)

    def activated(self, venv):
        """this environment with venv's activate script sourced.  An install's activate does more than a bare
        venv's: a source build's and a release's both put the install's lib/ on LD_LIBRARY_PATH, which some of its
        binaries (vcfwave) cannot run without"""
        with self._env_lock:
            if venv not in self._activated:
                # the script as $1, not $0: virtualenv's activate refuses to run if $0 is itself
                r = subprocess.run(['bash', '-c', '. "$1" > /dev/null && env -0', 'bigeval',
                                    os.path.join(venv, 'bin', 'activate')], capture_output=True, stdin=subprocess.DEVNULL)
                if r.returncode:
                    raise RuntimeError('could not source {}/bin/activate: {}'.format(venv, r.stderr.decode().strip()))
                env = dict(kv.split('=', 1) for kv in r.stdout.decode(errors='replace').split('\0') if '=' in kv)
                # bash's own bookkeeping, not activate's
                for k in ('_', 'SHLVL', 'PWD', 'OLDPWD'):
                    env.pop(k, None)
                    if k in os.environ:
                        env[k] = os.environ[k]
                self._activated[venv] = env
            return dict(self._activated[venv])

    def env(self, venv=None, cactus_bin=True):
        """the environment to run a step in: venv activated and, unless cactus_bin is False, this install's
        binaries next on the PATH.  With cactus_bin False, every cactus install is taken off the PATH, so a tool
        that needs its own bcftools (say) does not pick up cactus's.  Toil workers inherit the leader's
        environment, so this is theirs too."""
        venv = venv or self.venv
        env = self.activated(venv)
        env.pop('PYTHONHOME', None)
        path = [d for d in env.get('PATH', '').split(':') if d and d != os.path.join(venv, 'bin') and not (
            os.path.exists(os.path.join(d, 'cactus_consolidated')) or os.path.exists(os.path.join(d, 'cactus')))]
        env['PATH'] = ':'.join([os.path.join(venv, 'bin')] + ([self.bin] if cactus_bin else []) + path)
        env['VIRTUAL_ENV'] = venv
        return env

    def step(self, name, fn, *args):
        """run fn(*args) as step name, unless an earlier attempt finished it"""
        marker = self.path('steps', name + '.json')
        m = read_json(marker)
        if m and m['status'] == 'done':
            return
        start = time.time()
        write_json(marker, {'status': 'running', 'start': start})
        self.log('{}: started'.format(name))
        try:
            fn(*args)
        except BaseException as e:
            write_json(marker, {'status': 'failed', 'start': start, 'end': time.time(), 'error': str(e)})
            self.log('{}: FAILED: {}'.format(name, e))
            raise
        write_json(marker, {'status': 'done', 'start': start, 'end': time.time()})
        self.log('{}: done after {}'.format(name, hms(time.time() - start)))

    def parallel(self, *fns):
        """call each of fns at once and wait for all of them -- a failure does not stop the others, which may be
        hours into a cluster run -- then raise the first failure"""
        with ThreadPoolExecutor(len(fns)) as ex:
            futures = [ex.submit(fn) for fn in fns]
        errors = [f.exception() for f in futures if f.exception() is not None]
        if errors:
            raise errors[0]
        return [f.result() for f in futures]

    def sh(self, cmd, name, env=None):
        """run cmd here, appending its output to logs/name.out"""
        out = self.path('logs', name + '.out')
        self.log('{}: {}'.format(name, ' '.join(shlex.quote(c) for c in cmd)))
        with open(out, 'a') as log:
            log.write('### {} {}\n'.format(now(), ' '.join(shlex.quote(c) for c in cmd)))
            log.flush()
            r = subprocess.run(cmd, stdout=log, stderr=log, stdin=subprocess.DEVNULL, cwd=self.out,
                               env=env or self.env())
        if r.returncode:
            raise RuntimeError('{} exited with {}, see {}'.format(name, r.returncode, out))

    def toil_opts(self):
        opts = ['--batchSystem', 'slurm', '--doubleMem', 'true', '--retryCount', '5']
        if self.p.get('partition'):
            opts += ['--slurmPartition', self.p['partition']]
        return opts + shlex.split(self.p.get('toil_opts') or '')

    def toil(self, name, cmd, output=None, env=None, toil_opts=None):
        """run the Toil workflow cmd, its jobstore (js/name) inserted as its first argument.  It is restarted when an
        earlier attempt left its jobstore.  A workflow removes its jobstore when it finishes, so if output exists
        and the jobstore does not, it finished and is not run again."""
        js = self.path('js', name)
        if output and os.path.exists(output) and not os.path.exists(js):
            self.log('{}: {} is there and its jobstore is not: already finished'.format(name, output))
            return
        restart = ['--restart'] if os.path.exists(js) else []
        env = env or self.env()
        # so stop can find the right toil to kill it with
        with open(js + '.toil', 'w') as f:
            f.write(shutil.which('toil', path=env['PATH']) or 'toil')
        self.sh([cmd[0], js] + cmd[1:] + ['--logFile', self.path('logs', name + '.log')] +
                (self.toil_opts() if toil_opts is None else toil_opts) + restart, name, env)

    def sbatch(self, name, script, cpus=1, mem='4G', time='24:00:00', venv=None, cactus_bin=True):
        """run the bash script as a slurm job and wait for it, in the environment of env(venv, cactus_bin), from DIR"""
        job = self.path('jobs', name + '.sh')
        log = self.path('logs', name + '.out')
        venv = venv or self.venv
        with open(job, 'w') as f:
            f.write('#!/bin/bash\nset -eo pipefail\n')
            f.write('source {}\n'.format(shlex.quote(os.path.join(venv, 'bin', 'activate'))))
            if cactus_bin:
                f.write('export PATH={}:$PATH\n'.format(shlex.quote(self.bin)))
            f.write('cd {}\n'.format(shlex.quote(self.out)))
            f.write(textwrap.dedent(script).strip() + '\n')
        cmd = ['sbatch', '--wait', '--parsable', '-J', self.tag + '.' + name, '-c', str(cpus), '--mem', mem,
               '--time', time, '-o', log]
        partition = self.partition(time)
        if partition:
            cmd += ['-p', partition]
        cmd += shlex.split(self.p.get('sbatch_opts') or '') + [job]
        self.log('{}: sbatch {} ({} cores, {})'.format(name, job, cpus, mem))
        r = subprocess.run(cmd, stdout=subprocess.DEVNULL, stdin=subprocess.DEVNULL,
                           env=self.env(venv, cactus_bin))
        if r.returncode:
            raise RuntimeError('slurm job {} failed, see {}'.format(name, log))

    def partition(self, limit):
        """the partition for a job with time limit limit: --partition, else the one Toil would choose -- one without
        GRES, then of the lowest priority tier, then with the shortest time limit that holds the job.  The cluster's
        default partition may be too short for it.  None leaves it to slurm."""
        if self.p.get('partition'):
            return self.p['partition']
        with self._slurm_lock:
            if self._partitions is None:
                r = subprocess.run(['sinfo', '-a', '-h', '-o', '%P %G %l %p'], capture_output=True, text=True,
                                   stdin=subprocess.DEVNULL)
                self._partitions = [l.split() for l in r.stdout.splitlines() if len(l.split()) == 4]
        fits = [(gres != '(null)', int(tier) if tier.isdigit() else sys.maxsize, slurm_seconds(t), name.rstrip('*'))
                for name, gres, t, tier in self._partitions if slurm_seconds(t) >= slurm_seconds(limit)]
        return min(fits)[3] if fits else None

    @contextlib.contextmanager
    def cache_lock(self, key):
        os.makedirs(self.cache_path('locks'), exist_ok=True)
        with open(self.cache_path('locks', re.sub(r'[^\w.-]', '_', key) + '.lock'), 'w') as f:
            fcntl.flock(f, fcntl.LOCK_EX)
            yield

    def cached(self, key, fn):
        """the cache directory key, made by fn(directory) unless a run already made it.  Runs sharing the cache
        wait for each other rather than make it twice."""
        d = self.cache_path(key)
        if not os.path.exists(d + '.done'):
            with self.cache_lock(key):
                if not os.path.exists(d + '.done'):
                    shutil.rmtree(d, ignore_errors=True)
                    os.makedirs(d)
                    fn(d)
                    open(d + '.done', 'w').close()
        return d

    def fetch(self, url):
        """url downloaded into the cache once, checked against its Content-Length: its local path"""
        def download(d):
            dst = os.path.join(d, os.path.basename(urllib.parse.urlparse(url).path))
            self.log('downloading {}'.format(url))
            for attempt in range(5):
                try:
                    with urllib.request.urlopen(url, timeout=600) as resp, open(dst, 'wb') as f:
                        size = resp.headers.get('Content-Length')
                        shutil.copyfileobj(resp, f, 1 << 24)
                    if size is not None and os.path.getsize(dst) != int(size):
                        raise IOError('got {} of {} bytes'.format(os.path.getsize(dst), size))
                    break
                except OSError as e:
                    if attempt == 4:
                        raise RuntimeError('could not download {}: {}'.format(url, e))
                    self.log('download of {} failed ({}); retrying'.format(url, e))
                    time.sleep(60)
            with open(os.path.join(d, 'URL'), 'w') as f:
                f.write(url + '\n')
        if not re.match(r'(https?|ftp)://', url):
            return url
        key = 'download/' + hashlib.sha1(url.encode()).hexdigest()[:16]
        return os.path.join(self.cached(key, download), os.path.basename(urllib.parse.urlparse(url).path))

    def fetch_all(self, urls, parallel=4):
        with ThreadPoolExecutor(parallel) as ex:
            return list(ex.map(self.fetch, urls))

    def steps(self):
        out = collections.OrderedDict()
        if os.path.isdir(self.path('steps')):
            for f in os.listdir(self.path('steps')):
                if f.endswith('.json'):
                    out[f[:-5]] = read_json(self.path('steps', f))
        return collections.OrderedDict(sorted(out.items(), key=lambda kv: kv[1]['start']))

    def hours(self, *names):
        """wall-clock hours steps names took, all together"""
        st = [self.steps().get(n) for n in names]
        if not st or any(s is None or s['status'] != 'done' for s in st):
            return ''
        return round((max(s['end'] for s in st) - min(s['start'] for s in st)) / 3600, 2)

def summarize(r):
    """the run's row: who and what, then the suite's metrics"""
    metrics = r.suite.summarize(r)
    ident = r.p['identity']
    row = collections.OrderedDict()
    steps = r.steps()
    row['date'] = datetime.date.fromtimestamp(max(s['end'] for s in steps.values())).isoformat()
    row['label'] = r.p['label']
    row['cactus_commit'] = ident['cactus_commit'][:12]
    row['cactus_version'] = ident['cactus_version']
    row['dirty'] = 'yes' if ident['dirty'] or (ident['checkout_head'] and
                            ident['cactus_commit'][:12] != ident['checkout_head'][:12]) else 'no'
    row['config'] = ','.join(r.p['set']) or ('stock' if not r.p['config'] else os.path.basename(r.p['config']))
    row['cactus_opts'] = r.p.get('cactus_opts') or ''
    row['eval'] = r.suite.EVAL_VERSION
    row.update(metrics)
    row['run'] = r.out
    summary = collections.OrderedDict([('suite', r.p['suite']), ('row', row), ('params', r.p), ('steps', steps)])
    write_json(r.path('summary.json'), summary)
    with open(r.path('summary.tsv'), 'w') as f:
        f.write('\t'.join(row) + '\n' + '\t'.join(str(v) for v in row.values()) + '\n')
    return row

def print_rows(rows, out=sys.stdout):
    """rows (dicts) as aligned columns, one column per row: they are wide"""
    keys = []
    for row in rows:
        keys += [k for k in row if k not in keys]
    width = max(len(k) for k in keys)
    for k in keys:
        out.write('{:{}s}  {}\n'.format(k, width, '  '.join('{:>14}'.format(str(row.get(k, ''))) for row in rows)))

def drive(out):
    """run (or carry on with) the run in out, here, in the foreground"""
    r = Run(out)
    lock = open(r.path('driver.lock'), 'w')
    try:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
    except OSError:
        sys.exit('another driver is running {}'.format(r.out))
    write_json(r.path('driver.json'), {'pid': os.getpid(), 'pgid': os.getpgid(0), 'host': socket.gethostname(),
                                       'start': time.time()})
    r.log('{} run {} ({}): cactus {} from {}'.format(r.p['suite'], r.out, r.p['label'],
                                                     r.p['identity']['cactus_commit'], r.venv))
    try:
        r.suite.run(r)
        row = summarize(r)
    except Exception:
        traceback.print_exc()
        r.log('the run stopped; once whatever failed is fixed, carry on with: {} resume {}'.format(sys.argv[0], r.out))
        return 1
    print_rows([row])
    r.log('finished; {} record {} adds this row to results/{}.tsv'.format(sys.argv[0], r.out, r.p['suite']))
    return 0

def detach(out):
    """carry on with the run in out in the background, in a session of its own, logging to out/driver.log"""
    log = os.path.join(out, 'driver.log')
    with open(log, 'a') as f:
        p = subprocess.Popen([sys.executable, os.path.abspath(__file__), 'resume', out], stdout=f,
                             stderr=subprocess.STDOUT, stdin=subprocess.DEVNULL, start_new_session=True)
    print('running in the background as pid {}, logging to {}'.format(p.pid, log))
    print('  {} status {}'.format(sys.argv[0], out))
    return 0

def cmd_run(a):
    out = os.path.abspath(a.out)
    if os.path.exists(out):
        sys.exit('{} already exists: resume it, or choose another --out'.format(out))
    suite = load_suite(a.suite)
    venv = os.path.abspath(a.cactus)
    if not os.path.exists(os.path.join(venv, 'bin', 'activate')):
        sys.exit('--cactus {}: not a virtualenv'.format(venv))
    if not os.path.isdir(os.path.join(os.path.dirname(venv), 'bin')):
        sys.exit('--cactus {}: no bin/ beside it for the cactus binaries'.format(venv))
    if not shutil.which('sbatch'):
        sys.exit('no sbatch on the PATH: run this on a node that can submit to slurm')
    ident = identity(venv)
    if ident['dirty']:
        print('warning: {} has local changes'.format(ident['checkout']), file=sys.stderr)
    if ident['checkout_head'] and ident['cactus_commit'] != ident['checkout_head']:
        print('warning: the installed cactus is {} but {} is at {}: the installed code is what runs'.format(
            ident['cactus_commit'], ident['checkout'], ident['checkout_head']), file=sys.stderr)
    params = collections.OrderedDict([
        ('suite', a.suite), ('label', a.label or os.path.basename(out)), ('cactus', venv),
        ('config', os.path.abspath(a.config) if a.config else None), ('set', a.set or []),
        ('cactus_opts', a.cactus_opts), ('cache', os.path.abspath(a.cache or os.environ.get('BIGEVAL_CACHE') or
                                                                  os.path.join(os.path.dirname(out), 'cache'))),
        ('partition', a.partition), ('sbatch_opts', a.sbatch_opts), ('toil_opts', a.toil_opts),
        ('created', now()), ('identity', ident), ('harness', harness_identity())])
    os.makedirs(out)
    for d in ('steps', 'logs', 'js', 'jobs'):
        os.makedirs(os.path.join(out, d))
    try:
        edit_config(params['config'] or stock_config(venv), params['set'], os.path.join(out, 'config.xml'))
        r = Run(out, params)
        suite.check(r)
    except Exception as e:
        shutil.rmtree(out)
        sys.exit('{}: {}'.format(a.suite, e))
    write_json(os.path.join(out, 'run.json'), params)
    return detach(out) if a.detach else drive(out)

def cmd_resume(a):
    out = os.path.abspath(a.dir)
    Run(out)
    return detach(out) if a.detach else drive(out)

def cmd_status(a):
    r = Run(a.dir)
    print('{} run {}: {}, cactus {}{}'.format(r.p['suite'], r.out, r.p['label'], r.p['identity']['cactus_commit'][:12],
                                              ' (with ' + ','.join(r.p['set']) + ')' if r.p['set'] else ''))
    d = read_json(r.path('driver.json'))
    # a running driver holds the lock
    with open(r.path('driver.lock'), 'a') as lock:
        try:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
            running = False
        except OSError:
            running = True
    print('driver: running, pid {} on {}'.format(d['pid'], d['host']) if running else 'driver: not running')
    for name, s in r.steps().items():
        # a step left running by a driver that is gone was cut off
        status = 'stopped' if s['status'] == 'running' and not running else s['status']
        took = hms((s.get('end') or time.time()) - s['start']) if status != 'stopped' else ''
        print('  {:24s} {:8s} {:>10s}{}'.format(name, status, took, '  ' + s['error'] if s.get('error') else ''))
    if os.path.exists(r.path('summary.tsv')):
        print()
        print_rows([read_json(r.path('summary.json'))['row']])

def cmd_stop(a):
    r = Run(a.dir)
    d = read_json(r.path('driver.json'))
    if d and d['host'] != socket.gethostname():
        sys.exit('the driver runs on {}: stop it from there'.format(d['host']))
    # each Toil workflow is told to stop, so that it takes its own jobs down with it
    for f in os.listdir(r.path('js')):
        if f.endswith('.toil') and os.path.exists(r.path('js', f[:-5])):
            with open(r.path('js', f)) as t:
                subprocess.run([t.read().strip(), 'kill', r.path('js', f[:-5])])
    if d:
        try:
            os.killpg(d['pgid'], signal.SIGTERM)
        except ProcessLookupError:
            pass
    # and the jobs it gave sbatch, which are named for the run
    if shutil.which('squeue'):
        jobs = subprocess.run(['squeue', '-h', '-u', os.environ['USER'], '-o', '%i %j'], capture_output=True,
                              text=True).stdout.split('\n')
        ids = [j.split()[0] for j in jobs if len(j.split()) == 2 and j.split()[1].startswith(r.tag + '.')]
        if ids:
            subprocess.run(['scancel'] + ids)
    print('stopped; {} resume {} carries on'.format(sys.argv[0], r.out))

def cmd_summarize(a):
    print_rows([summarize(Run(a.dir))])

def cmd_record(a):
    r = Run(a.dir)
    summary = read_json(r.path('summary.json'))
    if summary is None:
        sys.exit('{} has no summary.json: has it finished?'.format(r.out))
    row = summary['row']
    table = a.results or os.path.join(RESULTS, r.p['suite'] + '.tsv')
    rows, header = [], []
    if os.path.exists(table):
        with open(table) as f:
            rd = csv.DictReader(f, delimiter='\t')
            header, rows = list(rd.fieldnames), list(rd)
    if any(x.get('run') == row['run'] and x.get('date') == row['date'] for x in rows):
        sys.exit('{} already has this run'.format(table))
    header += [k for k in row if k not in header]
    rows.append(row)
    os.makedirs(os.path.dirname(os.path.abspath(table)), exist_ok=True)
    with open(table + '.tmp', 'w') as f:
        w = csv.DictWriter(f, header, delimiter='\t', lineterminator='\n', restval='')
        w.writeheader()
        w.writerows(rows)
    os.replace(table + '.tmp', table)
    print('added {} to {}'.format(row['label'], table))

def cmd_compare(a):
    runs = [Run(d) for d in a.dirs]
    suites = {r.p['suite'] for r in runs}
    if len(suites) != 1:
        sys.exit('those runs are of different suites: {}'.format(', '.join(sorted(suites))))
    rows = [(read_json(r.path('summary.json')) or {}).get('row') or {'label': r.p['label']} for r in runs]
    print_rows(rows)
    if hasattr(runs[0].suite, 'compare'):
        print()
        runs[0].suite.compare(runs, sys.stdout)

def cmd_list(a):
    for name in suite_names():
        s = load_suite(name)
        print('{:12s} {}'.format(name, s.DESCRIPTION))

def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest='cmd', required=True)
    sub.add_parser('list', help='the suites')
    s = sub.add_parser('run', help='start a run')
    s.add_argument('suite')
    s.add_argument('--cactus', required=True, help='the cactus virtualenv to test (its binaries in the bin/ beside it)')
    s.add_argument('--out', required=True, help='the run directory, which must not exist yet')
    s.add_argument('--label', help='a name for the run in its row [the run directory\'s name]')
    s.add_argument('--config', help='the cactus config to start from [the install\'s own]')
    s.add_argument('--set', action='append', metavar='[ELEMENT.]ATTR=VALUE',
                   help='set a config attribute, which must be in the config already, on one element (repeatable)')
    s.add_argument('--cactusOpts', dest='cactus_opts', help='more options for the suite\'s cactus command, in quotes')
    s.add_argument('--cache', help='the cache shared between runs [$BIGEVAL_CACHE, else cache/ beside the run]')
    s.add_argument('--partition', help='slurm partition, for both Toil and sbatch [chosen by each job\'s time limit, '
                   'as Toil does]')
    s.add_argument('--sbatchOpts', dest='sbatch_opts', help='more options for every sbatch, in quotes')
    s.add_argument('--toilOpts', dest='toil_opts', help='more options for every Toil workflow, in quotes')
    s.add_argument('--detach', action='store_true', help='run in the background')
    s = sub.add_parser('resume', help='carry on with a run that stopped')
    s.add_argument('dir')
    s.add_argument('--detach', action='store_true', help='run in the background')
    for name, desc in (('status', 'where a run is'), ('stop', 'stop a run and its cluster jobs'),
                       ('summarize', 'make a finished run\'s row again')):
        sub.add_parser(name, help=desc).add_argument('dir')
    s = sub.add_parser('record', help='add a finished run\'s row to its suite\'s results table')
    s.add_argument('dir')
    s.add_argument('--results', help='the table [results/SUITE.tsv beside this script]')
    s = sub.add_parser('compare', help='finished runs of one suite side by side, in detail')
    s.add_argument('dirs', nargs='+')
    a = ap.parse_args()
    sys.exit(globals()['cmd_' + a.cmd](a) or 0)

if __name__ == '__main__':
    main()
