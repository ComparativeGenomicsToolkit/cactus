# Cactus
[![Build Status](https://travis-ci.org/ComparativeGenomicsToolkit/cactus.svg?branch=master)](https://travis-ci.org/ComparativeGenomicsToolkit/cactus)

Cactus is a reference-free whole-genome multiple alignment program. The principal algorithms are described here: https://doi.org/10.1101/gr.123356.111 

## Acknowledgements

Cactus uses many different algorithms and individual code contributions, principally from Joel Armstrong, Glenn Hickey, Mark Diekhans and Benedict Paten. We are particularly grateful to:

-  Yung H. Tsin and Nima Norouzi for contributing their 3-edge connected components program code, which is crucial in constructing the cactus graph structure, see: Tsin,Y.H., "A simple 3-edge-connected component algorithm," Theory of Computing Systems, vol.40, No.2, 2007, pp.125-142.
- Bob Harris for providing endless support for his LastZ pairwise, blast-like genome alignment tool.


## Setup
### System requirements
Cactus uses substantial resources. For primate-sized genomes (3 gigabases each), you should expect Cactus to use approximately 120 CPU-days of compute per genome, with about 120 GB of RAM used at peak. The requirements scale roughly quadratically, so aligning two 1-megabase bacterial genomes takes only 1.5 CPU-hours and 14 GB RAM.

Note that to run even the very small evolverMammals example, you will need 2 CPUs and 12 GB RAM. The actual resource requirements are much less, but the individual jobs have resource estimates based on much larger alignments, so the jobs will refuse to run unless there are enough resources to meet their estimates.
### Virtual environment
To avoid problems with conflicting versions of dependencies on your system, we strongly recommend installing Cactus inside a Python [virtual environment](https://virtualenv.pypa.io/en/stable/). Note that Cactus will only currently work with Python 2.7, until some of our dependencies become Python 3 compatible.

To install the `virtualenv` command, if you don't have it already, run:
```
pip install virtualenv
```

To set up a virtual environment in the directory `cactus_env`, run:
```
virtualenv cactus_env
```

Then, to enter the virtualenv, run:
```
source cactus_env/bin/activate
```

You can always exit out of the virtualenv by running `deactivate`. The rest of the README assumes you're running inside a virtual environment.
### Install Cactus and its dependencies
Cactus uses [Toil](http://toil.ucsc-cgl.org/) to coordinate its jobs. To install Toil into your environment, run:
```
pip install --upgrade toil[all]
```

Finally, to install Cactus, from the root of the `cactus` repository, run:
```
pip install --upgrade .
```

### Install Docker or Singularity
By default, Cactus uses Docker to run its compiled components (to avoid making you install dependencies). It can instead use Singularity to run its binaries, or use a locally installed copy. To select a different way of running the binaries, you can use the `--binariesMode singularity` or `--binariesMode local` options.

Note, you must use Singularity 2.3 - 2.6 or Singularity 3.1.0+. Singularity 3 versions below 3.1.0 are incompatible with cactus (see [issue #55](https://github.com/ComparativeGenomicsToolkit/cactus/issues/55) and [issue #60](https://github.com/ComparativeGenomicsToolkit/cactus/issues/60)).

### Compile Cactus executables (if not using Docker/Singularity)
By default Cactus uses containers to distribute its binaries, because compiling its dependencies can sometimes be a pain. If you can use Docker or Singularity, you can skip this section. However, in some environments (e.g. HPC clusters) you won't be able to use Docker or Singularity, so you will have to compile the binaries and install a few dependencies.

First, ensure you have KyotoTycoon installed. If you have root access, it is available through most package managers under `kyototycoon` or `kyoto-tycoon`. To compile it manually, you are best off using the [unofficial Altice Labs repository](https://github.com/alticelabs/kyoto). If you've installed KyotoTycoon (and its library, KyotoCabinet) from a package manager, you should be OK to go. If you've installed it in a non-standard location, however, (because you don't have root access, for example) you will need to set the following environment variables:
```
ttPrefix=<path of the PREFIX where you installed Kyoto>
export kyotoTycoonIncl="-I${ttPrefix}/include -DHAVE_KYOTO_TYCOON=1"
export kyotoTycoonLib="-L${ttPrefix}/lib -Wl,-rpath,${ttPrefix}/lib -lkyototycoon -lkyotocabinet -lz -lbz2 -lpthread -lm -lstdc++"
```
and copy the `ktserver` binary to somewhere on your PATH, and depending on your install directory, you may also need to add `${ttPrefix}/lib` to your LD_LIBRARY_PATH. (This can be a bit of a pain--we have an updated `scons`-based build system in the works that will automate most of this, but it's not ready yet.)

Once you have KyotoTycoon installed, you should be able to compile Cactus and its dependencies by running:
```
git submodule update --init
make
```

To run using these local executables, you will need to provide the `--binariesMode local` option to all `cactus` commands and add the `bin` directory to your PATH.
## System/cluster requirements
Cactus will take about 20 CPU-hours per bacterial-sized (~4 megabase) genome, about 20 CPU-days per nematode-sized (~100 megabase) genome, and about 120 CPU-days per mammalian-sized (~3 gigabase) genome. You will need at least one machine with very large amounts of RAM (150+ GB) to run mammalian-sized genomes. The requirements will vary a bit depending on how closely related your genomes are, so these are only rough estimates.
## Running
To run Cactus, the basic format is:
```
cactus <jobStorePath> <seqFile> <outputHal>
```

The `jobStorePath` is where intermediate files, as well as job metadata, will be stored. It must be accessible to all worker systems.

When first testing out Cactus on a new system or cluster, before running anything too large, try running the small (5 600kb genomes) simulated example in `examples/evolverMammals.txt`. It should take less than an hour to run on a modern 4-core system. That example, even though it's small, should be enough to expose any major problems Cactus may have with your setup. The command you should run is:
```
cactus jobStore examples/evolverMammals.txt examples/evolverMammals.hal --root mr
```
### Choosing how to run the Cactus binaries (Docker/Singularity/local)
By default, Cactus uses Docker to run its compiled components (to avoid making you install dependencies). It can instead use Singularity to run its binaries, or use a locally installed copy. To select a different way of running the binaries, you can use the `--binariesMode singularity` or `--binariesMode local` options. (If running using local binaries, you will need to make sure cactus's bin directory is in your `PATH`.)
### seqFile: the input file
The input file, called a "seqFile", is just a text file containing the locations of the input sequences as well as their phylogenetic tree. The tree will be used to progressively decompose the alignment by iteratively aligning sibling genomes to estimate their parents in a bottom-up fashion. Polytomies in the tree are allowed, though the amount of computation required for a sub-alignment rises quadratically with the degree of the polytomy. The file is formatted as follows:

    NEWICK tree (optional)
    name1 path1
    name2 path2
    ...
    nameN pathN

An optional * can be placed at the beginning of a name to specify that its assembly is of reference quality. This implies that it can be used as an outgroup for sub-alignments. If no genomes are marked in this way, all genomes are assumed to be of reference quality. The star should only be placed on the name-path lines and not inside the tree.

* The tree must be on a single line. All leaves must be labeled and these labels must be unique. Ancestors may be named, or left blank (in which case the ancestors in the final output will automatically be labeled Anc0, Anc1, etc.) Labels must not contain any spaces.
* Branch lengths that are not specified are assumed to be 1.
* Lines beginning with # are ignored. 
* Sequence paths must point to either a FASTA file or a directory containing 1 or more FASTA files.
* Sequence paths must not contain spaces.
* Each name / path pair must be on its own line
* `http://`, `s3://`, etc. URLs may be used.


Please ensure your genomes are *soft*-masked with RepeatMasker. We do some basic masking as a preprocessing step to ensure highly repetitive elements are masked when repeat libraries are incomplete, but genomes that aren't properly masked can still take tens of times longer to align that those that are masked. Hard-masking (totally replacing repeats with stretches of Ns) isn't necessary, and is strongly discouraged (you will miss a *lot* of alignments!).

Example:

	  # Sequence data for progressive alignment of 4 genomes
	  # human, chimp and gorilla are flagged as good assemblies.
	  # since orang isn't, it will not be used as an outgroup species.
     (((human:0.006,chimp:0.006667):0.0022,gorilla:0.008825):0.0096,orang:0.01831);
     *human /data/genomes/human/human.fa
     *chimp /data/genomes/chimp/
     *gorilla /data/genomes/gorilla/gorilla.fa
     orang /cluster/home/data/orang/
### Running locally
There isn't much to configure if running locally. Most importantly, if on a shared system, you can adjust the maximum number of processors used with `--maxCores <N>` (by default, Cactus will use all cores).
### Running on a cluster
Cactus (through Toil) supports many batch systems, including LSF, SLURM, GridEngine, Parasol, and Torque. To run on a cluster, simply add `--batchSystem <batchSystem>`, e.g. `--batchSystem gridEngine`. If your batch system needs additional configuration, Toil exposes some [environment variables](http://toil.readthedocs.io/en/3.10.1/developingWorkflows/batchSystem.html#batch-system-enivronmental-variables) that can help.
### Running on the cloud
Cactus supports running on AWS, Azure, and Google Cloud Platform using [Toil's autoscaling features](https://toil.readthedocs.io/en/latest/running/cloud/cloud.html). For more details on running in AWS, check out [these instructions](doc/running-in-aws.md) (other clouds are similar).
## Using the output
Cactus outputs its alignments in the [HAL](https://github.com/ComparativeGenomicsToolkit/hal) format. This format represents the alignment in a reference-free, indexed way, but isn't readable by many tools. To export a MAF (which by its nature is usually reference-based), you can use the `hal2maf` tool to export the alignment from any particular genome: `hal2maf <hal> --refGenome <reference> <maf>`.

You can use the alignment to generate gene annotatations for your assemblies, using the [Comparative Annotation Toolkit](https://github.com/ComparativeGenomicsToolkit/Comparative-Annotation-Toolkit).
# Frequently Asked Questions
Q: I'm running under macOS using the Docker functionality and get an error from Docker: `docker: Error response from daemon: Mounts denied: [...]`

A: Go to your Docker preferences. In the "File Sharing" tab, double-click the last entry ("/path/to/exported/directory") and type in `/var/folders`. (Don't use the `+` button, it won't work because it resolves symlinks before adding).

The reason you have to do this is that the Docker VM requires explicitly listing the directories that can be bind-mounted. The default temp directory on macOS (`/var/folders/...`) is *symlinked* to a directory that is already listed as bind-mountable, but Docker checks the listing before resolving the symlink, returning an error.
