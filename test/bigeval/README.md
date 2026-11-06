# Big evaluations

Whole-genome tests of cactus that are too big for CI: each aligns a fixed set of public genomes on a slurm
cluster with the cactus install under test, scores the result against outside truth, and boils the scores down to
one row of numbers. Run them by hand before a release or to judge a major change. The rows are kept here in
`results/`, one table per suite, with the cactus commit each came from.

Nothing in here is data: genomes, annotations and truth sets are downloaded from their public sources.

| suite | what | truth | time | disk |
|---|---|---|---|---|
| `busco5` | `cactus` on human, mouse, opossum, chicken and gar | BUSCO gene models (compleasm), ENCODE CTCF cCREs | about a day, mostly queueing | 10 GB |
| `hprc30` | `cactus-pangenome` on CHM13, GRCh38 and 28 HPRC haplotypes, as in HPRC v2.1 | HPRC dipcall calls, GIAB HG002 Q100 | 5-20 h to build, 1 h to score | 130 GB |

Each suite's module, `suites/<suite>.py`, says what every column of its row means, and how much of a difference is
noise. In short:

- `busco5`: CPU hours, then the recall and precision of BUSCO coding bases aligned in phase to their ortholog
  (`busco_recall`, `busco_precision`), overall and at each depth of the tree (gar, the root, first; then mouse,
  opossum, chicken); human CTCF sites aligned, and kept beyond chance; genome-wide coverage; `cactus_consolidated`'s
  time and peak memory.
- `hprc30`: CPU hours, then the recall, precision and F1 of the graph's raw VCF (sites over 10 kb dropped) against
  the 14 samples' dipcall calls, pooled; errors (FN + FP) of it and of the vcfwave VCF, against dipcall and against
  GIAB's HG002 small and structural variant benchmarks; graph size, path fragments, snarl sizes and the share of
  input left out of the graph.

## Running

On the cluster's head node: it has to be able to submit to slurm and reach the internet. The driver and the Toil
leaders run there, and everything else goes to slurm, so it needs neither AVX2 nor much memory.

```
B=/shared/work/bigeval    # anywhere on the cluster's shared filesystem
test/bigeval/bigeval.py run busco5 --cactus /path/to/cactus/venv-cactus --out $B/busco5-$(date +%F)-master --detach
test/bigeval/bigeval.py status $B/busco5-2026-10-21-master
```

`--cactus` is the virtualenv of the install to test, inside the install, the way both a source build and a binary
release lay it out: its binaries are in the `bin/` beside it. The installed `cactus` package's commit is what the
row records. A venv installed without `-e` keeps the code it was installed with, so the checkout can move on, but
don't rebuild or reinstall that cactus while a run is using it: jobs pick up binaries and code as they start.

Options worth knowing:

- `--set [ELEMENT.]ATTR=VALUE` (repeatable) changes the install's own config, e.g. `--set baseAligner=abpoa`. The
  attribute must already be there, on exactly one element (of tag `ELEMENT`, if given), so a typo or a build without
  the option fails at once, not days later. `--config FILE` starts from another config instead.
- `--cactusOpts '...'` adds options to the suite's cactus command.
- `--dropOpts=FLAGS` leaves flags (options without a value) out of the suite's cactus command, e.g.
  `--dropOpts=--mgSplit`. A flag the suite does not pass fails the run at once.
- `--label NAME` names the run in its row (default: the run directory's name).
- `--partition`, `--sbatchOpts`, `--toilOpts` for the cluster.
- `--cache DIR`: see below. Default `$BIGEVAL_CACHE`, else `cache/` beside the run directory.

A run that stopped (a dead head node, a failed job) carries on with `bigeval.py resume DIR [--detach]`: finished
steps are skipped and Toil workflows restarted from their jobstores. `bigeval.py stop DIR` stops the driver, its Toil
workflows (`toil kill`) and its slurm jobs.

## Results

A finished run writes its row to `DIR/summary.tsv`, and the row plus everything behind it (parameters, what was
tested, each step's times) to `DIR/summary.json`. To keep it:

```
test/bigeval/bigeval.py record DIR        # appends the row to test/bigeval/results/<suite>.tsv
git add test/bigeval/results/<suite>.tsv  # and commit
```

`record` can run anywhere that can read the run directory. When a suite gains columns, `summarize DIR` remakes a
finished run's row and `record --replace DIR` puts it in place of the old one. GitHub shows the tables as tables,
scrolling sideways. To see runs side by side in detail (per pair, per
ancestor, per sample, and so on, depending on the suite):

```
test/bigeval/bigeval.py compare DIR1 DIR2 ...
```

A row reads from the left: which run, then the suite's headline numbers (CPU hours, recall and precision), then
the suite's details, then exactly what was run:

| column | |
|---|---|
| `label` | `--label` |
| `cactus_version` | `git describe` of the install's checkout, or the installed commit if it is not a checkout |
| ... | the suite's columns, its headline numbers first; `suites/<suite>.py` says what each means |
| `config` | the `--set` changes, `stock` for none, or the `--config` file's name |
| `cactus_opts` | `--cactusOpts`, then `without FLAG` for each `--dropOpts` flag |
| `date` | the day the run started |
| `cactus_commit` | the installed cactus package's commit |
| `dirty` | `yes` if the checkout had local changes or was not at the installed commit |
| `eval` | the suite's `EVAL_VERSION`: rows scored by different versions are not comparable |
| `harness_commit` | the commit of this harness the run was started from |
| `run` | the run directory |

## The cache

Runs share a cache of whatever does not depend on the install under test: the downloads, the BUSCO calls, the CTCF
motif scans and the scoring tools. Runs sharing a cache wait for each other rather than make anything twice. Keep
it: besides saving time it keeps rows comparable, because compleasm does not make the same BUSCO calls twice
(`busco_truth` in busco5's rows says which calls a row was scored against).

## Layout

```
bigeval.py      the driver: running, resuming, summarizing and recording runs
suites/         a module per suite, and the seqfiles and truth manifests it ships
scoring/        scorers the suites run on the cluster
results/        a table per suite, one row per recorded run
```

A run directory:

```
run.json        parameters, and what was tested
config.xml      the cactus config it aligned with
steps/          each step's status and times
logs/           a log per step: S.log is Toil's --logFile, S.out everything else
js/  jobs/      Toil jobstores, and the scripts given to sbatch
summary.tsv  summary.json
```

## Adding a suite

A suite is a module in `suites/` with `NAME`, `DESCRIPTION`, `EVAL_VERSION`, and three functions of the run
(`bigeval.Run`):

- `check(r)`: raise if the install or config cannot run the suite (called before anything is launched)
- `run(r)`: the pipeline, as steps: `r.step(name, fn, ...)` runs `fn` unless an earlier attempt finished it, and
  `r.parallel(...)` runs steps at once. `r.toil` runs a Toil workflow (restarting it on resume), `r.sbatch` a job,
  `r.fetch` a download into the cache and `r.cached` anything else into the cache.
- `summarize(r)`: the metrics, an ordered dict

and optionally `compare(runs, out)`. Bump `EVAL_VERSION` whenever a change to the suite or its scorers changes the
numbers.
