"""primates5: progressive cactus on five catarrhines, human to macaque, scored as busco5 is.

busco5's ancestors are 0.33 to 0.73 apart (in the scaled tree cactus measures divergence in), so it says nothing about
the shallow ancestors that make up most of a big alignment: 72% of the 447-way's are under 0.05.  These five are a node
apart, from human-chimp (0.013) to the catarrhine root (0.085), so settings that depend on divergence, like <pecan
shallowRepeatMaskDivergence>, can be set on either side of each node with --set.  Five primates also align in a
fraction of busco5's time.

The scoring is busco5's, run on its own copy of busco5's code with the genomes below; see busco5.py for what every
column means and how much of a difference is noise.  The headline numbers are cpu_h, busco_recall and busco_precision,
then macaque_recall and macaque_precision, the same at the root.  The pairs of each split are scored on their own:
chimp (human-chimp), gorilla (the Homininae node), orangutan (the Hominidae node) and macaque (the root).  Between
genomes this close BUSCO recall and precision are near their ceilings, so read nc_inlocus and ctcf_beyond_chance too:
they are what sees non-coding sequence, repeats included, aligned outside its ortholog.
"""
import collections, importlib.util, os

NAME = 'primates5'
DESCRIPTION = 'progressive cactus on human, chimp, gorilla, orangutan and macaque: BUSCO, CTCF and coverage'
# bump when anything below or in busco5's scoring changes the numbers, so rows scored differently are not compared
EVAL_VERSION = 1

# busco5's code, loaded as a module of our own so that setting its constants here leaves busco5 itself as it is
_spec = importlib.util.spec_from_file_location('bigeval_suite_primates5_busco5',
                                               os.path.join(os.path.dirname(os.path.abspath(__file__)), 'busco5.py'))
_busco5 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_busco5)

_busco5.NAME = NAME
_busco5.EVAL_VERSION = EVAL_VERSION
_busco5.SEQFILE = 'primates5.seqfile'
_busco5.CHROMINFO = 'primates5.chrom-info'
# as in busco5: the 577-way's (ROADIES) branch lengths run shorter than cactus expects, so they are doubled
CACTUS_OPTS = _busco5.CACTUS_OPTS
# human is GRCh38, as in busco5, so the CTCF cCREs (and the BUSCO calls and motif scans in the cache) are shared
_busco5.CTCF_SOURCE = 'GCA_000001405.15'
# a pair is filed under the deepest split it spans, as in busco5: human-chimp, then each further genome with
# everything closer to human, so gorilla's pairs show the Homininae node, orangutan's the Hominidae node and macaque's
# the root
_busco5.DEPTH = collections.OrderedDict([('GCA_000001405.15', 0), ('GCF_028858775.2', 0), ('GCF_029281585.2', 1),
                                         ('GCF_028885655.2', 2), ('GCF_037993035.1', 3)])
_busco5.GROUPS = ['chimp', 'gorilla', 'orangutan', 'macaque']
_busco5.COMMON = {'GCA_000001405.15': 'human', 'GCF_028858775.2': 'chimp', 'GCF_029281585.2': 'gorilla',
                  'GCF_028885655.2': 'orangutan', 'GCF_037993035.1': 'macaque'}

check = _busco5.check
run = _busco5.run
summarize = _busco5.summarize
compare = _busco5.compare
