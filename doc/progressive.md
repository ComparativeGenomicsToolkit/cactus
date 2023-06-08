# Progressive Cactus

Progressive Cactus is included in the [Cactus Software Package](../README.md).  It can be used to align hundreds of vertebrate genomes.  A phylogenetic tree is required as input in order to define the progressive decomposition. For aligning samples from the *same* species, use the [Minigraph-Cactus](./pangenome.md) pipeline instead. 

Please cite the [Progressive Cactus paper](https://doi.org/10.1038/s41586-020-2871-y) when using Cactus.

## Table of Contents

* [Quick-start](#quick-start)
* [Interface](#interface)
* [Using the HAL Output](#using-the-hal-output)
* [Using Docker](#using-docker)
* [Running in the Cloud](#running-on-the-cloud)
* [Running on a Cluster](#running-on-a-cluster)
* [Running Step-by-step](#running-step-by-step)
* [Running on Terra or Cromwell](#running-on-terra-or-cromwell)
* [Updating Alignments](#updating-alignments)
* [GPU Acceleration](#gpu-acceleration)
* [Pre-Alignment Checklist](#pre-alignment-checklist)
* [Frequently Asked Questions](#frequently-asked-questions)

## Quick Start

A small test example is included in `cactus/examples/evolverMammals.txt`.  To align these genomes, install the [latest release](https://github.com/ComparativeGenomicsToolkit/cactus/releases) and run from the cactus directory:

```
cactus ./js ./examples/evolverMammals.txt ./evolverMammals.hal
```

The alignment can be inspected with
```
halStats ./evolverMammals.hal

hal v2.2
((simHuman_chr6:0.144018,((simMouse_chr6:0.084509,simRat_chr6:0.091589)mr:0.171974,simGorilla:0.075)AncGorilla:0.1)Anc1:0.020593,(simCow_chr6:0.18908,simDog_chr6:0.16303)Anc2:0.032898)Anc0;

GenomeName, NumChildren, Length, NumSequences, NumTopSegments, NumBottomSegments
Anc0, 2, 535128, 13, 0, 17165
Anc1, 2, 561672, 7, 19963, 24668
simHuman_chr6, 0, 601863, 1, 25791, 0
AncGorilla, 2, 578003, 4, 26318, 57158
mr, 2, 607607, 5, 55240, 60749
simMouse_chr6, 0, 636262, 1, 62021, 0
simRat_chr6, 0, 647215, 1, 61761, 0
simGorilla, 0, 599081, 1, 58887, 0
Anc2, 2, 573075, 19, 20396, 61873
simCow_chr6, 0, 602619, 1, 61912, 0
simDog_chr6, 0, 593897, 1, 61666, 0
```

or converted to MAF with
```
cactus-hal2maf ./js evolverMammals.hal evolverMammals.maf.gz --refGenome simHuman_chr6 --chunkSize 1000000
```

## Interface

*Note*: See the [step-by-step interface](#running-step-by-step) to see how to run Progressive Cactus one job at a time.

Cactus takes as input a set of **softmasked** genome assemblies in fasta format (they can be gzipped) and a phylogenetic tree relating them. Each input genome should correspond to a leaf in the tree. It outputs a multiple genome alignment in [HAL](https://github.com/ComparativeGenomicsToolkit/hal), which includes ancestral sequences, one for each internal node in the tree.   

To run Cactus, the basic format is:
```
cactus <jobStorePath> <seqFile> <outputHal>
```

The `jobStorePath` is where intermediate files, as well as job metadata, [will be kept by Toil](https://toil.readthedocs.io/en/latest/running/introduction.html#job-store). **It must be accessible to all worker systems.** 

The `seqFile` is a text file containing the locations of the input sequences as well as their phylogenetic tree. The tree will be used to progressively decompose the alignment by iteratively aligning sibling genomes to estimate their parents in a bottom-up fashion. Polytomies in the tree are allowed, though the amount of computation required for a sub-alignment rises quadratically with the degree of the polytomy.  Cactus uses the predicted branch lengths from the tree to determine appropriate pairwise alignment parameters, allowing closely related species to be aligned more quickly with no loss in accuracy. The file is formatted as follows:

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
* FASTA files may be gzipped (since Cactus v2.2.0)
* Sequence paths must not contain spaces.
* Each name / path pair must be on its own line
* `http://`, `s3://`, etc. URLs may be used.

Please ensure your genomes are *soft*-masked with RepeatMasker. We do some basic masking as a preprocessing step to ensure highly repetitive elements are masked when repeat libraries are incomplete, but genomes that aren't properly masked can still take tens of times longer to align that those that are masked. Hard-masking (totally replacing repeats with stretches of Ns) isn't necessary, and is strongly discouraged (you will miss a *lot* of alignments!).

An example seqfile can be found [here](../examples/evolverMammals.txt).

## Using the HAL Output

The `outputHal` file represents the multiple alignment, including all input and inferred ancestral sequences.  It is stored in HAL format, and can be accessed with [HAL tools](https://github.com/ComparativeGenomicsToolkit/Hal), which are all included in Cactus either as static binaries for the binary release, or within the Docker image for the Docker release.

Please [cite HAL](https://doi.org/10.1093/bioinformatics/btt128).

### MAF

The HAL format represents the alignment in a reference-free, indexed way, but isn't readable by many tools. Many applications will require direct access to alignment columns, which requires "transposing" the row-based HAL format. This is best done by converting to [MAF](https://genome.cse.ucsc.edu/FAQ/FAQformat.html#format5), which can be computationally expensive for larger files. A toil-powered distributed MAF exporter is provided with Cactus.  Note that MAF is a reference-based format, so a reference (which can be any genome, leaf or ancestor, in the HAL file must be specified).

```
cactus-hal2maf ./js ./evolverMammals.hal evolverMammals.maf.gz --refGenome simHuman_chr6 --chunkSize 1000000 --noAncestors
```

`cactus-hal2maf`, in addition to providing parallelism with Toil, also adds [TAFFY](https://github.com/ComparativeGenomicsToolkit/taffy)-based normalization and is therefore recommended over running `hal2maf` directly. To better understand the normalzation options, see the [TAFFY](https://github.com/ComparativeGenomicsToolkit/taffy) link above. `cactus-hal2maf` provides an interface to toggle them, and attempts to use sensible defaults. They can be completely deactivated with the `--raw` option. 

The various batching options can be used to tune distributed runs on very large inputs. For example, to run 4 batches, each on a 32-core EC2 node but only processing 8 chunks with `taffy add-gap-bases` in parallel, these options could be used
```
--chunkSize 1000000 --batchCount 4 --batchCores 32 --batchParallelTaf 8 --batchSystem mesos --provisioner aws --defaultPreemptable --nodeStorage 2000 --maxNodes 4 --nodeTypes r5.8xlarge 
```

Depending on the application, you may want to handle duplication events differently when creating the MAF. Three different modes are available via the `--dupeMode` option.

* "single" : Uses greedy heuristics to pick the copy for each species that results in fewest mutations and block breaks. Recommended when visualizing via BigMaf (see below)
* "ancestral" : Restricts the duplication relationships shown to only those orthologous to the reference genome according to the HAL tree. There may be multiple orthologs per genome. This relies on the dating of the duplication in the hal tree (ie in which genome it is explicitly self-aligned) and is still a work in progress. For example, in a tree with `((human,chimp),gorilla)`, if a duplication in human is collapsed (ie a single copy) in the human-chimp ancestor, then it would not show up on the human-referenced MAF using this option. But if the duplication is not collapsed in this ancestor (presumably because each copy has an ortholog in chimp and gorilla), then it will be in the MAF because the duplication event was higher in the tree.
* "all" : (default) All duplications are written, including ancestral events (orthologs) and paralogs in the reference. 

Usually a reference genome is specified with `--refGenome` and ancestral genomes are excluded `--noAncestors`. Since the default reference is the root of the alignment, `--noAncestors` can only be specified if a leaf genome is used with `--refGenome`. 

### BigMaf

[UCSC BigMaf](https://genome.ucsc.edu/goldenPath/help/bigMaf.html) is an indexed version of MAF (see above) that can viewed on the Genome Browser. `cactus-maf2bigmaf` can be used to convert MAF (as output by `cactus-hal2maf`) into BigMaf.

It is recommended to create the maf using the `--noAncestors --dupeMode single` options with `cactus-hal2maf`.

```
cactus-maf2bigmaf ./js ./evolverMammals.maf.gz ./evolverMammals.bigmaf.bb --refGenome simHuman_chr6 --halFile evolverMammals.hal
```

This will produce the BigMaf file `evolverMammals.bigmaf.bb` along with the summary file `evolverMammals.bigmaf.summary.bb` which is used by the browser for zoomed out summary display.

The chromosome sizes of the reference genome must be provided as input either directly with `--chromSizes` or via the original hal file via `--halFile`.

### Chains

The [UCSC Chain Format](https://genome.ucsc.edu/goldenPath/help/chain.html) is a concise way to represent pairwise alignments, and is used by the Genome Browser and some of its tools. HAL files can be converted into sets of Chain files using `cactus-hal2chain`.

For example
```
cactus-hal2chain ./js ./evolverMammals.hal chains-dir --refGenome simHuman_chr6 
```

will create `./chains-dir` and populate it with a Chain alignment between simHuman and each other leaf genome in evolverMammals.hal.

By default, chains will be created using `halLiftover` [as in CAT](https://github.com/ComparativeGenomicsToolkit/Comparative-Annotation-Toolkit/blob/fc1623da5df1309d2e2f0b9bb0363aaab84708f4/cat/chaining.py#L96-L98). An option `--useHalSynteny` is provided to use that tool instead.

See here for an all-vs-all script to make chains, including BigChain conversion: https://github.com/human-pangenomics/HPRC_Assembly_Hub/blob/main/chains/wdl/snakesonachain.wdl

### CAT

You can use the alignment to generate gene annotatations for your assemblies, using the [Comparative Annotation Toolkit](https://github.com/ComparativeGenomicsToolkit/Comparative-Annotation-Toolkit).

### Assembly Hubs

HAL alignments can also be displayed in the UCSC genome browser via creation of assembly hubs as described [here](https://github.com/ComparativeGenomicsToolkit/hal#displaying-in-the-ucsc-genome-browser-using-assembly-hubs).  `hal2assemblyHub.py` is included in Cactus. 

### Phast

Conservation scores can be computed using [phast](http://compgen.cshl.edu/phast/) either directly from the HAL (`halPhyloP`) or from the MAF. The phast binaries are included in the Cactus releases. 

## Using Docker

The Cactus Docker image contains everything you need to run Cactus (python environment, all binaries, system dependencies). For example, to run the test data:

```
docker run -v $(pwd):/data --rm -it quay.io/comparative-genomics-toolkit/cactus:v2.5.2 cactus /data/jobStore /data/evolverMammals.txt /data/evolverMammals.hal
```

Or you can proceed interactively by running
```
docker run -v $(pwd):/data --rm -it quay.io/comparative-genomics-toolkit/cactus:v2.5.2 bash
cactus /data/jobStore /data/evolverMammals.txt /data/evolverMammals.hal

```

## Running on the cloud

Cactus supports running on AWS using [Toil's autoscaling features](https://toil.readthedocs.io/en/latest/running/cloud/cloud.html). For more details on running in AWS, check out [these instructions](./running-in-aws.md).

Cactus can also be run on Google Cloud Platform via [Terra](#running-on-cromwell/terra).

## Running on a cluster

Cactus (through Toil) supports many batch systems in theory, including LSF, SLURM, GridEngine, Parasol, and Torque. To run on a cluster, add `--batchSystem <batchSystem>`, e.g. `--batchSystem gridEngine`. If your batch system needs additional configuration, Toil exposes some [environment variables](http://toil.readthedocs.io/en/3.10.1/developingWorkflows/batchSystem.html#batch-system-enivronmental-variables) that can help.

IMPORTANT:  It is highly recommend that one **not** run Cactus in its default mode using the Toil Grid Engine-like batch systems (GridEngine, HTCondor, LSF, SLURM, or Torque).  Cactus creates a very large number of small jobs, which can overwhelm these systems.  The work-around described [here](#running-the-step-by-step-workflow-direclty-in-toil) for clusters with large compute nodes available must be used instead.  **Update**:  Cactus version >= 2.0 with GPU enabled will spawn far fewer jobs which, in theory, should make this less of an issue.

## Running step by step

Breaking Cactus up into smaller jobs can be practical, both for development and debugging, and managing larger workflows.  Here is an example of how to break the Evolver Mammals example up into three steps: 1) Preprocessing 2) Blast 3) Multiple Aligment:
```
cactus-prepare examples/evolverMammals.txt --outDir steps-output --outSeqFile steps-output/evovlerMammals.txt --outHal steps-output/evolverMammals.hal --jobStore jobstore
```

It will print the sequence of commands to run the alignment step-by-step.  Blocks of commands within each alignment run can be run in parallel

`cactus-prepare` can also be used to simplify preprocessing sequences without decomposing the remaining workflow:

```
cactus-prepare examples/evolverMammals.txt --outDir steps-output --outSeqFile steps-output/evovlerMammals.txt --outHal steps-output/evolverMammals.hal --jobStore jobstore --preprocessOnly
```

`cactus-prepare-toil` shares the interface of `cactus-prepare` except instead of printing the command lines or WDL script, it runs them directly from Toil.  An example use case of this is within UCSC's kubernetes cluster.  Like many computing environments, the number of jobs that can be scheduled is limited, so running `cactus` directly using Toil's `kubernetes` batch system will swamp the cluster.  But if the computation can be broken up into a handful of steps, and a job is only created for each step (as in the Cromwell/WDL method), then it can run through.  So `cactus-prepare-toil` will run as a high-level Toil workflow on the specified batch system, and it will launch jobs for each command (`cactus-preprocess, cactus-blast, cactus-align`), and each one of these jobs will get scheduled on a node and run its command with the `singleMachine` batch system.  Here is an example invocation for kubernetes:

```
cactus-prepare-toil aws:us-west-2:<JOBSTORE-NAME> examples/evolverMammals.txt --binariesMode singularity --batchSystem kubernetes --outHal s3://<BUCKET-NAME>/out.hal --defaultDisk 20G --defaultMemory 12G --defaultCores 4
```

## Running on Terra or Cromwell 

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

Here is an example of some settings that have worked on a mammalian-sized genome alignment on Terra.  It's important to align the [resources](https://cromwell.readthedocs.io/en/stable/RuntimeAttributes) requested (CPU and memory) to N1 instance types as found [here](https://gcloud-compute.com/instances.html).  Note that disk will be rounded up to the [nearest multiple of 375G](https://cromwell.readthedocs.io/en/stable/RuntimeAttributes/#disks), with only [these multiples](https://github.com/ComparativeGenomicsToolkit/cactus/issues/792) supported: `[0, 1, 2, 3, 4, 5, 6, 7, 8, 16, 24]`. In the example below this is:

* **preprocess/blast**: `n1-standard-32 + 8 v100 GPUs`
* **align**: `n1-highmem-64` (can lower to `n1-highmem-32` using `--alignCores 32 --alignMemory 208Gi` for shorter branches/smaller genomes)
* **append**: `n1-highmem-16`

```
cactus-prepare --wdl mammals.txt --noLocalInputs --preprocessBatchSize 5 \
               --preprocessDisk 375Gi --preprocessCores 32 --preprocessMemory 120Gi \
               --blastDisk 375Gi --blastCores 32 --gpu 8 --blastMemory 120Gi \
               --alignDisk 375Gi --alignCores 64 --alignMemory 416Gi \
               --halAppendDisk 3000Gi  --defaultMemory 104Gi  > mammals.wdl
```

If the workflow fails for whatever reason, it can be edited (to, say, increase job requirements) then resumed as follows:
* In the Workflows tab, click the scripts link beside "Source:" to go back to the Firecloud page to edit the WDL script
* Edit it and "Save a New Snapshot"
* Back in the Terra Workflows tab for the workflow, refresh the page, and select the new snapshot from the "Snapshots" menu.
* Click the "Save" button, ensure that "Use call caching" is ticked, then "Run Analysis" again to resume the workflow.

**Important note on resuming Terra workflows**

The above instructions to use Terra's call caching do not work reliably anymore.  This is really frustrating, as you can be many days and dollars into a workflow, need to adjust the WDL for whatever reason, and Terra will ignore all intermediate files and restart from scratch when you rerun for reasons only it understands. But if the results are in your GCP bucket somewhere (which they will be as long as you are not explicitly removing them), you can still use them by editing the WDL to incorporate them.  This can be done with `cactus-terra-helper`. You can use the Terra interface (by clicking on pretty much any file) to find the root bucket prefix of all intermediate files of a given run.  Once you have that, run

```
gsutil ls -r gs://<BUCKET/PREFIX> | cactus-terra-helper resume mammals.wdl > mammals-resume.wdl
```

The output script will remove all WDL calls for which output files were found in the bucket, and replace references to them to full paths of the intermediate files.

The same script can be used to download all the logs off Terra, which can be useful.  This command (the `-l` is important) will download the latest version of each log to the present directory. 

```
gsutil ls -l -r gs://<BUCKET/PREFIX> | cactus-terra-helper scrape-logs
```

## Updating Alignments

Cactus supports incrementally updating existing alignments to add, remove, or update genomes. The process involves minor surgery on the output HAL files. See [this document](updating-alignments.md) for details. [cactus-update-prepare](cactus-update-prepare.md) can be used to simplify this process!

## GPU Acceleration

[SegAlign](https://github.com/ComparativeGenomicsToolkit/SegAlign), a GPU-accelerated version of lastz, can be used in the "preprocess" and "blast" phases to speed up the runtime considerably, provided the right hardware is available. Unlike lastz, the input sequences do not need to be chunked before running SegAlign, so it also reduces the number of Toil jobs substantially.  The [GPU-enabled Docker releases](https://github.com/ComparativeGenomicsToolkit/cactus/releases) have SegAlign turned on by default and require no extra options from the user.  Otherwise, it is possible to [manually install it](https://github.com/ComparativeGenomicsToolkit/SegAlign#-dependencies) and then enable it in `cactus` using the `--gpu` command line option. One effective way of ensuring that only GPU-enabled parts of the workflow are run on GPU nodes is on Terra with `cactus-prepare --gpu --wdl` (see above example).

By default `--gpu` will give all available GPUs to each SegAlign job. This can be tuned by passing in a numeric value, ex `--gpu 8` to assign 8 GPUs to each SegAlign job.  In non-single-machine batch systems, it is mandatory to set an exact value with `--gpu`.  

GPUs must
* support CUDA
* have at least 8GB GPU memory (for mammal-sized input)
* have 1-2 CPU cores available each for spawning lastz jobs

We've tested SegAlign on Nvidia V100 and A10G GPUs. See the Terra example above for suggested node type on GCP.   

Please [cite SegAlign](https://doi.ieeecomputersociety.org/10.1109/SC41405.2020.00043).


## Pre-Alignment Checklist

* Are the input sequences softmasked with RepeatMasker? For mammals we expect at least 40% of the genome to be masked this way. 
* Have you run a small [test alignment](../examples/evolverMammals.txt) to make sure Cactus is properly installed?
* Do you have at least one outgroup species?

## Frequently Asked Questions

**Q**: I'm running under macOS using the Docker functionality and get an error from Docker: `docker: Error response from daemon: Mounts denied: [...]`

**A**: Go to your Docker preferences. In the "File Sharing" tab, double-click the last entry ("/path/to/exported/directory") and type in `/var/folders`. (Don't use the `+` button, it won't work because it resolves symlinks before adding).

The reason you have to do this is that the Docker VM requires explicitly listing the directories that can be bind-mounted. The default temp directory on macOS (`/var/folders/...`) is *symlinked* to a directory that is already listed as bind-mountable, but Docker checks the listing before resolving the symlink, returning an error.

**Q**: Why does cactus write some files to `/tmp` when I specify another location with `--workDir`.  

**A**: This is a bug, but you can work around it by running `export TMPDIR=<your desired path>` before running cactus. 

**Q**: How exactly does Cactus use the branch lengths in the input tree?

**A**: The branch lengths are expected to be denoted in substitutions per site and are used in three ways. 

1) to determine lastz parameters. The pairwise distance between species is measured using the branch lengths, and mapped to a set of lastz parameters using the `<divergences>` (inside `<constants>`) and `<divergence>` (inside `<blast>`) elements in the [configuration XML](../src/cactus/cactus_progressive_config.xml). Faster parameters are used for more closely-related species. If you are aligning human and chimp with the correct branch lengths, their distance will be about 0.02, and it will use the fasest parameters which will be several times faster than if, say, the default branch length of 1 was used. 

2) to estimate ancestral bases. Here the relative branch lengths are more important -- the base of a descendant that is much nearer to the ancestor will provide more information and will be wieghted higher when estimating it.  

3) to calculate outgroups.  Branch lengths are taken into account by the greedy heuristic used to find the nearest outgroup to the given ancestral event. 

If you do not know the branch lengths, you can leave them out and Cactus will use its default (1). This will cause the alignment to be slower than necessary but the results shouldn't be affected much.  Otherwise, even inexact branch lengths should be fine.  For closely related species `mash dist` is a very easy way to estimate a pairwise distance (`mash` is now included in Cactus). Often you can use `mash` and already-published trees to come up with your branch lengths. 

We are currently working on incorporating a fast genome tree estimation workflow with Cactus.  

**Q**: How much do I have to worry about uncertainty in my guide tree?

**A**:  Changes in the tree topology will affect the alignment, but like the Progressive Cactus paper shows, it'll be fairly minimal for small, local changes. For most data, there is no single correct tree anyway, due to things like incomplete lineage sorting, hybridization and lateral gene transfer.

That said, cactus can handle multifurcations up to a point: runtime increases quadratically with the number of species in the split. So depending on your genome size, you probably won't be able to go much beyond 5. This should be more accurate in theory, but I haven't run the experiments to back it up.

**Q**: I'm running out of memory, or getting crashes, or very long runtimes in one of the `paf_xxxx` tools (`paf_tile, paf_to_bed` etc.). What can I do?

**A**: This is almost always due to Cactus having found too many pairwise alignments in the all-to-all lastz mapping (blast) phase. The only way to get around this is my softmasking the input genomes before running Cactus. For most species, this is best done with RepeatMasker. We do intend to work on lifting this requirement in the future by making cactus's own repeatmasking more robust. 

**Q**: The `--gpu` option isn't working for me.

**A**: Unless you've set up SegAlign yourself, the GPU option will only work using the gpu-enabled Docker image (name ending in `-gpu`).  If you are running directly from the container make sure you use `docker run`'s `--gpus` option to enable GPUs in your container. If you are using `singularity`, the option is `--nv`.

**Q**: But what if I want to use `--gpu` on my cluster? When I try from inside the GPU-enabled container, none of my cluster commands (ex `qsub`) are available.

**A**: Install the Cactus binary release as described in the instructions on the Releases page.  But run Cactus with `--binariesMode docker` (or `singularity`).  This will let Cactus run SegAlign (and all other binaries) directly from the container, while itself running from the Python virtualenv.

**Q**: I get an error to the effect of `ERROR: No matching distribution found for toil[aws]==xxxx` when trying to install Toil.

**A**: This is probably happening because you are using Python 3.6. Toil and Cactus require Python >= 3.7.  Use `python3 --version` to check your Python version.
