# Cactus
[![Build Status](https://travis-ci.org/ComparativeGenomicsToolkit/cactus.svg?branch=master)](https://travis-ci.org/ComparativeGenomicsToolkit/cactus)

Cactus is a reference-free whole-genome multiple alignment program. The principal algorithms are described here: https://doi.org/10.1101/gr.123356.111 

Please subscribe to the [cactus-announce](https://groups.google.com/d/forum/cactus-announce) low-volume mailing list so we may reach about releases and other announcements.

## Acknowledgements

Cactus uses many different algorithms and individual code contributions, principally from Joel Armstrong, Glenn Hickey, Mark Diekhans and Benedict Paten. We are particularly grateful to:

- Yung H. Tsin and Nima Norouzi for contributing their 3-edge connected components program code, which is crucial in constructing the cactus graph structure, see: Tsin,Y.H., "A simple 3-edge-connected component algorithm," Theory of Computing Systems, vol.40, No.2, 2007, pp.125-142.
- Bob Harris for providing endless support for his [LastZ](https://github.com/lastz/lastz) pairwise, blast-like genome alignment tool.
- Sneha Goenka and Yatish Turakhia for the [GPU-accelerated version of LastZ](https://github.com/ComparativeGenomicsToolkit/SegAlign).


## Setup

### System requirements
We regularly test on Ubuntu 18.04 (Bionic) and to a more limited degree on Mac OS X (using Docker).

Cactus requires Python 3.

Cactus uses substantial resources. For primate-sized genomes (3 gigabases each), you should expect Cactus to use approximately 120 CPU-days of compute per genome, with about 120 GB of RAM used at peak. The requirements scale roughly quadratically, so aligning two 1-megabase bacterial genomes takes only 1.5 CPU-hours and 14 GB RAM.

Note that to run even the very small evolverMammals example, you will need 2 CPUs and 12 GB RAM. The actual resource requirements are much less, but the individual jobs have resource estimates based on much larger alignments, so the jobs will refuse to run unless there are enough resources to meet their estimates.

IMPORTANT:  It is highly recommend that one **not** run Cactus using the Toil Grid Engine-like batch systems (GridEngine, HTCondor, LSF, SLURM, or Torque).  Cactus creates a very large number of small jobs, which can overwhelm these systems.

### Installation Overview

There are many different ways to install and run Cactus:
* [Docker Image](#docker-image)
* [Precompiled Binaries](#precompiled-binaries)
* [Build From Source](#build-from-source)
* [Python Install with Docker Binaries](#python-install-with-docker-binaries)

#### Docker Image

Cactus docker images are hosted on [quay](https://quay.io/repository/comparative-genomics-toolkit/cactus).  The image for the latest release is listed on the [Releases Page](https://github.com/ComparativeGenomicsToolkit/cactus/releases).  Here is an command line to run the included evolver mammals example with release 1.1.1
```
wget https://raw.githubusercontent.com/ComparativeGenomicsToolkit/cactus/master/examples/evolverMammals.txt
docker run -v $(pwd):/data --rm -it quay.io/comparative-genomics-toolkit/cactus:v1.1.1 cactus /data/jobStore /data/evolverMammals.txt /data/evolverMammals.hal --root mr --binariesMode local

```

#### Precompiled Binaries

Precompiled binaries can be found on the [Releases Page](https://github.com/ComparativeGenomicsToolkit/cactus/releases).  Download by clicking the `cactus-bin-vx.x.x.tar.gz` and install following the instructions in the `BIN-INSTALL.md` link.

#### Build From Source

##### Install the Cactus Python package and its dependencies

To avoid problems with conflicting versions of dependencies on your system, we strongly recommend installing Cactus inside a Python 3 [virtual environment](https://virtualenv.pypa.io/en/stable/). To install the `virtualenv` command, if you don't have it already, run:
```
python3 -m pip install virtualenv
```

To set up a virtual environment in the directory `cactus_env`, run:
```
python3 -m virtualenv -p python3.6 cactus_env
```

Then, to enter the virtualenv, run:
```
source cactus_env/bin/activate
```

You can always exit out of the virtualenv by running `deactivate`.


To install Cactus in Python, clone it and **its submodules with --recursive** from github and install it with pip:
```
git clone https://github.com/ComparativeGenomicsToolkit/cactus.git --recursive
cd cactus
pip install --upgrade setuptools pip
pip install --upgrade -r toil-requirement.txt
pip install --upgrade .
```

##### Build the Cactus Binaries

Several binaries are required to run Cactus.  They can be built as follows:

Compile time settings can be overridden by creating a make include file in the top level cactus directory.  
```
cactus/include.local.mk
```

Cactus has several dependencies that need to be installed on the system, including HDF5. HDF5 is available through most package managers (`apt-get install libhdf5-dev`) or can be manual installed from source files at [The HDF Group](https://www.hdfgroup.org/).   HDF5 should be configured with the `--enable-cxx` option. If you've installed it in a non-standard location, have the `h5c++` command in your `PATH` or add this to `include.local.mk`:
```
export PATH := <hdf5 bin dir>:${PATH}
```

You can use the the [Dockerfile](Dockerfile) as a guide to see how all dependencies are installed with `apt` on Ubuntu.

In the top level cactus directory.  The binaries can then be built with
```
make -j $(nproc)
```
and added to the PATH with
```
export PATH=$(pwd)/bin:$PATH
```

To use HAL python scripts such as `hal2mafMP.py`, add the submodules directory to the PYTHONPATH with
```
export PYTHONPATH=$(pwd)/submodules:$PYTHONPATH
```

#### Python Install With Docker Binaries

Cactus can be setup and used in a virtual environment as in the [previous section](#build-from-source), without compiling the binaries.  When used like this (which will happen automatically when running `cactus` without the appropriate binaries in the `PATH` environment variable), a Docker image will be automatically pulled to run commands as needed.  The main use case for this is running with Toils AWS provisioner as [described here](doc/running-in-aws.md).

Singularity binaries can be used in place of docker binaries with the `--binariesMode singularity` flag.  Note, you must use Singularity 2.3 - 2.6 or Singularity 3.1.0+. Singularity 3 versions below 3.1.0 are incompatible with cactus (see [issue #55](https://github.com/ComparativeGenomicsToolkit/cactus/issues/55) and [issue #60](https://github.com/ComparativeGenomicsToolkit/cactus/issues/60)).

The `--binariesMode local` flag can be used to force `cactus` to run local binaries -- this is the default behavior if they are found. 

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

Within an hour at most (on modern computers), you should have a [HAL](https://github.com/ComparativeGenomicsToolkit/hal) file which relates simulated mouse and rat genomes.


### seqFile: the input file
The input file, called a "seqFile", is just a text file containing the locations of the input sequences as well as their phylogenetic tree. The tree will be used to progressively decompose the alignment by iteratively aligning sibling genomes to estimate their parents in a bottom-up fashion. Polytomies in the tree are allowed, though the amount of computation required for a sub-alignment rises quadratically with the degree of the polytomy.  Cactus uses the predicted branch lengths from the tree to determine appropriate pairwise alignment parameters, allowing closely related species to be aligned more quickly with no loss in accuracy. The file is formatted as follows:

    NEWICK tree
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

### Running on a cluster
Cactus (through Toil) supports many batch systems, including LSF, SLURM, GridEngine, Parasol, and Torque. To run on a cluster, simply add `--batchSystem <batchSystem>`, e.g. `--batchSystem gridEngine`. If your batch system needs additional configuration, Toil exposes some [environment variables](http://toil.readthedocs.io/en/3.10.1/developingWorkflows/batchSystem.html#batch-system-enivronmental-variables) that can help.

### Running on the cloud
Cactus supports running on AWS, Azure, and Google Cloud Platform using [Toil's autoscaling features](https://toil.readthedocs.io/en/latest/running/cloud/cloud.html). For more details on running in AWS, check out [these instructions](doc/running-in-aws.md) (other clouds are similar).

### Running step by step

#### Printing a list of commands to run locally

Breaking Cactus up into smaller jobs can be practical, both for development and debugging, and managing larger workflows.  Here is an example of how to break the Evolver Mammals example up into three steps: 1) Preprocessing 2) Blast 3) Multiple Aligment:
```
cactus-prepare examples/evolverMammals.txt --outDir steps-output --outSeqFile steps-output/evovlerMammals.txt --outHal steps-output/evolverMammals.hal --jobStore jobstore
```

It will print the sequence of commands to run the alignment step-by-step.  Blocks of commands within each alignment run can be run in parallel

`cactus-prepare` can also be used to simplify preprocessing sequences without decomposing the remaining workflow:

```
cactus-prepare examples/evolverMammals.txt --outDir steps-output --outSeqFile steps-output/evovlerMammals.txt --outHal steps-output/evolverMammals.hal --jobStore jobstore --preprocessOnly
```

#### Creading a WDL script to run on Cromwell or Terra

The `--wdl` option in `cactus-prepare` can be used to generate a bespoke [WDL](https://github.com/openwdl/wdl/blob/master/versions/1.0/SPEC.md) script for running the alignment from the input seqFile.  Here is an example on how to run locally in [Cromwell](https://github.com/broadinstitute/cromwell)
```
cactus-prepare examples/evolverMammals.txt --wdl > evolver.wdl
wget https://github.com/broadinstitute/cromwell/releases/download/49/cromwell-49.jar
javac -jar ./cromwell-49.jar run evolver.wdl

```

To run on [Terra](https://terra.bio/), use the `--noLocalInputs` option to make sure no local files are embedded in the script.  Also, care must be taken to specify some minimum resource requirements.

```
cactus-prepare examples/evolverMammals.txt --wdl --noLocalInputs --alignCores 2 --defaultMemory 16G > evolver_terra.wdl

```

Then in Terra's [workspace menu](https://app.terra.bio/#workspaces):
* Create a new workspace if necessary with the "+" button
* Click on the workspace
* Click on the "DATA" tab in the workspace menu and use the "Files" link to upload `examples/evolverMammals.txt` to Goggle Cloud
* Click on the "WORKFLOWS" tab
* Click the "+" button to add a workflow
* Click the link in the bottom right to the "Broad Methods Repository"
* Click the "Create New Method... +" button
* Choose and namespace and name, then either upload or paste `evolver_terra.wdl` as created above and click "Upload"
* If this WDL is valid, you can use the "Export To Workspace" button to link it to the Terra Workspace (using a blank configuration)
* You can select the option to go back to the Terra Workspace, otherwise the workflow should now appear as a card in the Terra "workflows" tab the next time you navigate there or refresh
* To run it, click the workflow then click the "INPUTS" tab, and select the `evolverMammals.txt` file in the Attribute field for Task=`cactus_prepare` Variable=`prep_seq_file`
* Tick "Run workflow with inputs defined by file paths"
* Save and click "RUN ANALYSIS"

In the evolver example, all input sequences are specified in public URLs.  If sequences are not specified as URLs in the seqfile, then they must be uploaded in similar fashion to how the evolverMammals.txt was uploaded and selected in the example above.

Here is an example of some settings that have worked on a mammalian-sized genome alignment on Terra:

```
cactus-prepare --wdl mammals.txt --noLocalInputs --preprocessBatchSize 5 --alignDisk 3000G --halAppendDisk 3000G --preprocessDisk 3000G --defaultDisk 1000G --defaultCores 64 --gpu --gpuCount 8 --defaultMemory 385G > mammals.wdl

```

If the workflow fails for whatever reason, it can be edited (to, say, increase job requirements) then resumed as follows:
* In the Workflows tab, click the scripts link beside "Source:" to go back to the Firecloud page to edit the WDL script
* Edit it and "Save a New Snapshot"
* Back in the Terra Workflows tab for the workflow, refresh the page, and select the new snapshot from the "Snapshots" menu.
* Click the "Save" button, ensure that "Use call caching" is ticked, then "Run Analysis" again to resume the workflow.

#### Running the step-by-step workflow direclty in Toil

`cactus-prepare-toil` shares the interface of `cactus-prepare` except instead of printing the command lines or WDL script, it runs them directly from Toil.  An example use case of this is within UCSC's kubernetes cluster.  Like many computing environments, the number of jobs that can be scheduled is limited, so running `cactus` directly using Toil's `kubernetes` batch system will swamp the cluster.  But if the computation can be broken up into a handful of steps, and a job is only created for each step (as in the Cromwell/WDL method), then it can run through.  So `cactus-prepare-toil` will run as a high-level Toil workflow on the specified batch system, and it will launch jobs for each command (`cactus-preprocess, cactus-blast, cactus-align`), and each one of these jobs will get scheduled on a node and run its command with the `singleMachine` batch system.  Here is an example invocation for kubernetes:

```
cactus-prepare-toil aws:us-west-2:<JOBSTORE-NAME> examples/evolverMammals.txt --binariesMode singularity --batchSystem kubernetes --realTimeLogging --outHal s3://<BUCKET-NAME>/out.hal --defaultDisk 20G --defaultMemory 12G --defaultCores 4
```

## GPU Acceleration

A [GPU-accelerated version of lastz](https://github.com/ComparativeGenomicsToolkit/SegAlign) can be used in the `blast` phase to speed up the runtime considerably, provided the right hardware is available. The easiest way to use it is on Terra with `cactus-prepare --gpu --wdl` (see above example).  The [GPU-enabled Docker releases](https://github.com/ComparativeGenomicsToolkit/cactus/releases) have this turned on by default.  It is also possible to [manually install it](https://github.com/ComparativeGenomicsToolkit/SegAlign#-dependencies) from git and then enable it in `cactus` via the `--configFile` cactus option.  A template can be found [here](https://github.com/ComparativeGenomicsToolkit/cactus/blob/master/src/cactus/cactus_progressive_config.xml), and modified to activate the GPU by setting `gpuLastz="true"` and `realign="0"` in the `<caf>` section. 

## Using the output
Cactus outputs its alignments in the [HAL](https://github.com/ComparativeGenomicsToolkit/hal) format. This format represents the alignment in a reference-free, indexed way, but isn't readable by many tools. To export a MAF (which by its nature is usually reference-based), you can use the `hal2maf` tool to export the alignment from any particular genome: `hal2maf <hal> --refGenome <reference> <maf>`.

You can use the alignment to generate gene annotatations for your assemblies, using the [Comparative Annotation Toolkit](https://github.com/ComparativeGenomicsToolkit/Comparative-Annotation-Toolkit).

You can also [convert the HAL alignment into a Pangenome Graph](https://github.com/ComparativeGenomicsToolkit/hal#pangenome-graph-export-gfa-and-vg).  `hal2vg` is now included in the Cactus Docker images and binary release. 

## Updating existing alignments
Cactus supports incrementally updating existing alignments to add, remove, or update genomes. The process involves minor surgery on the output HAL files. See [this document](doc/updating-alignments.md) for details.
# Frequently Asked Questions
Q: I'm running under macOS using the Docker functionality and get an error from Docker: `docker: Error response from daemon: Mounts denied: [...]`

A: Go to your Docker preferences. In the "File Sharing" tab, double-click the last entry ("/path/to/exported/directory") and type in `/var/folders`. (Don't use the `+` button, it won't work because it resolves symlinks before adding).

The reason you have to do this is that the Docker VM requires explicitly listing the directories that can be bind-mounted. The default temp directory on macOS (`/var/folders/...`) is *symlinked* to a directory that is already listed as bind-mountable, but Docker checks the listing before resolving the symlink, returning an error.
