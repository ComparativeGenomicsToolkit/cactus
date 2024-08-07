# Minigraph-Cactus Pangenome Construction and Downstream Applications

## Table of Contents

* [Abstract](#abstract)
* [Key Reference Material](#key-reference-material)
* [Part 1: Pangenome Graph Construction](#part-1-pangenome-graph-construction)
     * [Cactus Setup](#cactus-setup)
     * [Input Data](#input-data)
     * [Build and Index the Pangenome Graph](#build-and-index-the-pangenome-graph)
     * [Running on Slurm](#running-on-slurm)     
* [Part 2: Pangenome Graph Properties](#part-2-pangenome-graph-properties)
     * [Basic Statistics](#basic-statistics)
     * [Subgraph Extraction](#subgraph-extraction)
     * [Format Conversion](#format-conversion)
     * [Panacus](#panacus)
* [Part 3: Mapping Reads to the Graph](#part-3-mapping-reads-to-the-graph)
     * [Short Read Mapping](#short-read-mapping)
     * [Long Read Mapping (stretch goal)](#long-read-mapping)
     * [Surjecting To BAM](#surjecting-to-bam)
* [Part 4: Genotyping and Variant Calling](#part-4-genotyping-and-variant-calling)
     * [Variant Calling with DeepVariant](#variant-calling-with-deepvariant)
     * [SV Genotyping with vg](#sv-genotyping-with-vg)
     * [SV Genotyping with pangenie (stretch goal)](#sv-genotyping-with-pangenie)

## Abstract

This is a tutorial written to support the **Reference Graph Pangenome Data Analysis Hackathon 2023 Nov. 13-17 in Cape Town, South Africa**. The aim is to provide detailed instructions on how to create a pangenome reference graph with Minigraph-Cactus then use it for some downstream analysis like variant calling and genotyping.

Unlike some previous workshops, and most of the existing Cactus documentation, this tutorial will focus on whole-genome human data. As such, it will need to be run over a period of time longer than a typical workshop session. The running times and memory usage of each command will be given wherever possible. 

**Update:** Many thanks to the workshop attendees for their patience and interest, and all credit to your feedback for helping to debug it! I will try to continue to support it at least for a little while.  Feel free to reach out to me on github or otherwise (even if you weren't in the workshop) with any questions or problems.

## Key Reference Material

Please visit these links for related material and background information before proceeding further. **The first link is essential and should absolutely be consulted before continuing, and the rest are highly recommended.** 

* [Minigraph-Cactus Manual](https://github.com/ComparativeGenomicsToolkit/cactus/blob/master/doc/pangenome.md): This is essential to read, and includes several small examples (with data) that should be run before tackling whole-genomes.
* [Minigraph-Cactus Paper](https://doi.org/10.1038/s41587-023-01793-w): The methods are described in detail here.
* [HPRC v1.1 Minigraph-Cactus Instructions](https://github.com/ComparativeGenomicsToolkit/cactus/blob/master/doc/mc-pangenomes/hprc-v1.1-mc.md): Commands and explanations in order to exactly reproduce the latest released HPRC graphs. The commands themselves assume a Slurm cluster but can be trivially modified to run on a single computer (remove `--batchSystem slurm`).
* [HPRC Graph Downloads](https://github.com/human-pangenomics/hpp_pangenome_resources/): Get the HPRC graphs here.
* [HPRC Paper](https://doi.org/10.1038/s41586-023-05896-x): Detailed analysis of the HPRC graph, and examples of many downstream applications of Minigraph-Cactus pangenomes. 
* @jeizenga's [2023 Memphis Workshop](https://github.com/pangenome/MemPanG23/blob/main/lessons/Day_3a_vg_mapping_and_calling.md), which served as an inspiration for this tutorial.

## Part 1: Pangenome Graph Construction

### Cactus Setup

**Important:** We will be using [Cactus v2.6.13](https://github.com/ComparativeGenomicsToolkit/cactus/releases/tag/v2.6.13) for this tutorial. Be warned that some steps may not work for older (or newer) versions.

For simplicity, all cactus will be run in "single-machine" mode via its [docker](https://www.docker.com/) image.  Cactus also supports distributed computing environments via [slurm](https://github.com/ComparativeGenomicsToolkit/cactus/blob/master/doc/progressive.md#running-on-a-cluster) and [AWS/Mesos](https://github.com/ComparativeGenomicsToolkit/cactus/blob/master/doc/running-in-aws.md). See the [Running on Slurm](#running-on-slurm) section for more details about running on a cluster.

In order to make sure `singularity` is working, try running the following and verify that you do not get an error. If this step does not work, you will need to consult your local sysadmin. 
```
singularity run docker://hello-world
```

### Input Data

As you've seen in the [Minigraph-Cactus Manual](https://github.com/ComparativeGenomicsToolkit/cactus/blob/master/doc/pangenome.md) (please go back and read it if you haven't already), the input is a list of sample name and genome assembly pairs.  For diploid assemblies, the convention of `SAMPLE.1 / SAMPLE.2` must be used (and dots avoided in sample names otherwise).

In addition to your samples of interest, you should include at least one reference genome. This will allow you to use reference coordinates to, for example, project variants on.  In this example, which is based on a small subset of 4 samples of the [HPRC]((https://github.com/ComparativeGenomicsToolkit/cactus/blob/master/doc/mc-pangenomes/hprc-v1.1-mc.md)) data, we will use GRCh38 and CHM13.

Please copy-paste the following data into `hprc10.seqfile`

```
GRCh38	https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/references/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta.gz
CHM13	https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/analysis_set/chm13v2.0.fa.gz
HG00438.1	https://s3-us-west-2.amazonaws.com/human-pangenomics/working/HPRC/HG00438/assemblies/year1_f1_assembly_v2_genbank/HG00438.paternal.f1_assembly_v2_genbank.fa.gz
HG00438.2	https://s3-us-west-2.amazonaws.com/human-pangenomics/working/HPRC/HG00438/assemblies/year1_f1_assembly_v2_genbank/HG00438.maternal.f1_assembly_v2_genbank.fa.gz
HG00621.1	https://s3-us-west-2.amazonaws.com/human-pangenomics/working/HPRC/HG00621/assemblies/year1_f1_assembly_v2_genbank/HG00621.paternal.f1_assembly_v2_genbank.fa.gz
HG00621.2	https://s3-us-west-2.amazonaws.com/human-pangenomics/working/HPRC/HG00621/assemblies/year1_f1_assembly_v2_genbank/HG00621.maternal.f1_assembly_v2_genbank.fa.gz
HG00673.1	https://s3-us-west-2.amazonaws.com/human-pangenomics/working/HPRC/HG00673/assemblies/year1_f1_assembly_v2_genbank/HG00673.paternal.f1_assembly_v2_genbank.fa.gz
HG00673.2	https://s3-us-west-2.amazonaws.com/human-pangenomics/working/HPRC/HG00673/assemblies/year1_f1_assembly_v2_genbank/HG00673.maternal.f1_assembly_v2_genbank.fa.gz
HG00733.1	https://s3-us-west-2.amazonaws.com/human-pangenomics/working/HPRC_PLUS/HG00733/assemblies/year1_f1_assembly_v2_genbank/HG00733.paternal.f1_assembly_v2_genbank.fa.gz
HG00733.2	https://s3-us-west-2.amazonaws.com/human-pangenomics/working/HPRC_PLUS/HG00733/assemblies/year1_f1_assembly_v2_genbank/HG00733.maternal.f1_assembly_v2_genbank.fa.gz

```

**If you are making a pangenome graph with your own data, this input listing should be the only part you need to change, but do see the explanation of the options below as some may require adjustments for different data sizes**. Also, nothing changes if you want to use haploid assemblies -- just do not use the `.1` and `.2` suffixes (see `CHM13` and `GRCh38` above).

### Build and Index the Pangenome Graph

I am going to run on 32-cores in order to simulate my understanding of an "average" node on your cluster. As you've seen in the [Minigraph-Cactus Manual](https://github.com/ComparativeGenomicsToolkit/cactus/blob/master/doc/pangenome.md) (please go back and read it if you haven't already), the simplest way to build the graph is with the `cactus-pangenome` command.

Here it is, with an explanation of each option following below.

**Update:** Previous versions of this document used `--giraffe clip filter` below.  Since Cactus v2.8.5, this is [no longer necessary](https://github.com/ComparativeGenomicsToolkit/cactus/pull/1424): you can use just `--giraffe filter` (or leave the `--giraffe` option out altogether to go all in on haplotype subsampling). Not using `--giraffe clip` drastically reduces peak memory usage.

```
rm -rf cactus-scratch && mkdir cactus-scratch

singularity exec -H $(pwd) docker://quay.io/comparative-genomics-toolkit/cactus:v2.6.13 \
cactus-pangenome ./js ./hprc10.seqfile --outDir ./hprc10 --outName hprc10 --reference GRCh38 CHM13 \
--filter 2 --haplo --giraffe filter --viz --odgi --chrom-vg clip filter --chrom-og --gbz clip filter full \
--gfa clip full --vcf --vcfReference GRCh38 CHM13 --logFile ./hprc10.log --workDir ./cactus-scratch \
--consCores 8 --mgMemory 128Gi
```

For `singularity exec`:
* `-H`: Set the home/working directory to the current directory
* `docker://quay.io/comparative-genomics-toolkit/cactus:v2.6.13`: the cactus docker image. It will be cached locally (probably in `~/.singularity` as a `.sif` file)

For `cactus-pangenome`:
* `./js`: Scratch directory that will be created for Toil's jobstore
* `./hprc10.seqfile`: The input samples and assemblies. This file was created above.
* `--outDir ./hprc10`: The output directory. All results will be here. 
* `--outName hprc10`: This will be the prefix of all output files.
* `--reference GRCh38 CHM13`: Specify these two samples as reference genomes. Reference samples are indexed a little different in `vg` to make their coordinates easier to use. Also, the first reference given (GRCh38 in this case), is used to anchor the entire graph and is treated differently than the other samples.  Please see [here for more details](https://github.com/ComparativeGenomicsToolkit/cactus/blob/master/doc/pangenome.md#reference-sample).
* `--filter 2`: Create an Allele-Frequency filtered (AF) graph that contains only nodes and edges supported by at least 2 haplotypes. This can lead to better mapping performance. `--filter 9` was used for the 90-assembly HPRC graph.  
* `--haplo`: We are actually phasing out the Allele-Frequency filtering as described above in favour or dynamic creation of personal pangenomes. Using this option will create the necessary indexes for this functionality.
* `--giraffe filter`: Make giraffe indexes for the Allele-Frequency filtered graph.
* `--viz`: Make an ODGI 1D visualization image for each chromosome.
* `--odgi`: Make an ODGI formatted whole-genome graph
* `--chrom-vg clip filter`: Make VG formatted chromosome graphs for the both the AF filtered and (default) clipped pangenome.
* `--chrom-og`: Make ODGI formatted chromosome graphs for the full (unclipped) graph. Useful for visualization. 
* `--gbz clip filter full`: Make GBZ formatted whole-genome graphs for AF filtered, (default) clipped and full (containing unaligned centromeres) graphs.
* `--gfa clip full`: Make GFA formatted whole-genome graphs for (default) clipped and full graphs.
* `--vcf`: Make a VCF (based on the first reference) version of the graph
* `--vcfReference GRCh38 CHM13`: Specify that we want two VCFs, one for each reference
* `--logFile ./hprc10.log`: All logging information will end up here in addition to `stderr`.  Important to save!
* `--consCores 8`: Specify 8 threads for each core cactus job (`cactus_consolidated`). By default it will use all cores available on your system.  By reducing to `8`, we attempt to run up to 4 chromosomes at once to save time (assuming 32 cores total). Note that this will increase peak memory usage.
* `--mgMemory 128Gi`: Override Cactus's estimated memory limit for `minigraph` construction to make sure it does not exceed what's available.  By default, Cactus is very conservative (estimates too much memory) in order to prevent jobs from being evicted from slurm clusters.  But we know here `126Gi` is fine. 
* `--workDir cactus-scratch` Location for cactus's temporary files.

All of the above is explained in more detail in the [Minigraph-Cactus Manual](https://github.com/ComparativeGenomicsToolkit/cactus/blob/master/doc/pangenome.md). We are erring on the side of producing lots of different indexes, but it's usually easier than going back and regenerating any forgotten ones.

Here are some details about the resources used. I'm on a big shared server using `docker run --cpus 32 --memory 250000000000` to emulate a smaller computer. This data is taken from `hpr10.log` which lists the wall time and memory usage of each command in the pipeline. The log for the full 90-way HPRC graph can be found [here](https://s3-us-west-2.amazonaws.com/human-pangenomics/pangenomes/freeze/freeze1/minigraph-cactus/hprc-v1.1-mc-grch38/hprc-v1.1-mc-grch38.log).

1) Minigraph Construction : 106Gi, ~3 hours
2) Minigraph Mapping : ~200Gi (max per-job 40Gi) ~2 hour
3) Cactus Alignment : ~65Gi (max per-job 16Gi) ~1 hour
4) Normalization and Indexing : ~64 Gi ~3 hours 
5) Overall : 11 hours (in environment with 32 cores / 256 Gb RAM)

Here are the output files:
```
4.0K	chrom-alignments
4.0K	chrom-subproblems
278M	hprc10.CHM13.raw.vcf.gz
1.6M	hprc10.CHM13.raw.vcf.gz.tbi
224M	hprc10.CHM13.vcf.gz
1.6M	hprc10.CHM13.vcf.gz.tbi
4.0K	hprc10.chroms
773M	hprc10.d2.dist
1.8G	hprc10.d2.gbz
35G	hprc10.d2.min
65M	hprc10.d2.snarls
1.1G	hprc10.dist
2.2G	hprc10.full.gbz
1.5G	hprc10.full.gfa.gz
9.7G	hprc10.full.hal
9.1G	hprc10.full.og
57M	hprc10.full.snarls
104M	hprc10.gaf.gz
1.8G	hprc10.gbz
769M	hprc10.gfa.fa.gz
1.4G	hprc10.gfa.gz
1.1G	hprc10.hapl
35G	hprc10.min
450M	hprc10.paf
7.3M	hprc10.paf.filter.log
117M	hprc10.paf.unfiltered.gz
279M	hprc10.raw.vcf.gz
1.6M	hprc10.raw.vcf.gz.tbi
907M	hprc10.ri
1.7K	hprc10.seqfile
52M	hprc10.snarls
406K	hprc10.stats.tgz
729M	hprc10.sv.gfa.gz
226M	hprc10.vcf.gz
1.6M	hprc10.vcf.gz.tbi
4.0K	hprc10.viz
```

There are four versions of the graph produced (please see [here](https://github.com/ComparativeGenomicsToolkit/cactus/blob/master/doc/pangenome.md#clipping-filtering-and-indexing) for more details), denoted by these preffixes:

* `hprc.sv` : This is the output of `minigraph` and contains only structural variants. The input haplotypes are not embedded as paths
* `hprc.full` : This is a base-level graph containing all sequence that could be assigned to a reference chromosome. Centromeres are included but are unaligned.
* `hprc.` : This is a subgraph of `hprc.full` but with centromeres removed.  This is usually the most relevant graph for analysis.
* `hprc.d2` : This is a subgraph of `hprc.` but with nodes and edges supported by fewer than 2 haplotypes removed. This graph yields better results for read mapping with the original `giraffe` pipeline. We used allele frequency filtering in the HPRC paper (`.d9` via ``-filter 9` for a `10%` cutoff) but have recently changed `giraffe` so that it is no longer necessary (more details later in the mapping section).

The graphs themselves are present in `.gfa.gz` (standard, text-based), `.gbz` (highly compressed, `vg`) and `.og` (`odgi`) formats. `vg giraffe` mapping requires the `.gbz` and `.hapl` index (or `.gbz`, `.dist` and `.min` for the original pipeline).

There are four VCF files, two each for GRCh38 and CHM13:

* `hprc10.raw.vcf.gz` and `hprc10.CHM13.raw.vcf.gz` : Output of `vg deconstruct` for the GRCh38- and CHM13-based graphs, respectively. These VCFs contain nested variants, and need special attention when using.
* `hprc10.vcf.gz` and `hprc10.CHM13.vcf.gz` : "Flattened" versions of the above VCFs (using `vcfbub`) that do not contain nested variants and will be more useful for standard tools.

For the HPRC, we took some extra normalization steps using `vcfwave` to realign the variants.  This gives a slightly cleaner VCF in some regions.  See [here](https://github.com/ComparativeGenomicsToolkit/cactus/blob/hprc-v1.1/doc/mc-pangenomes/hprc-v1.1-mc.md#vcf-postprocessing) for details.

There are four directories:

* `chrom-subproblems` : The per-chromosome *inputs* to `cactus`.  This directory has some useful statistics about the chromosome decomposition such as `contig_sizes.tsv` which shows the amount of sequence from each sample for each reference chromosome as well as `minigraph.split.log` which lists which chromosome each contig gets assigned to and why, along with all contigs excluded from further analysis (counted `_AMBIGUOUS_`) because they didn't align anywhere.  
* `chrom-alignments` : The raw, per-chromosome output of `cactus` including `.vg` and `.hal` files.  Only useful for debugging and / or re-running the last part of the pipeline (indexing and normalization).
* `hprc10.chroms` : Chromosome graphs in `.vg` and `.og` format.  Useful for debugging and visualization. If you are using `GRCh38` as a reference, the unplaced contigs will all get lumped into the `chrOther` graph.  
* `hprc10.viz` : ODGI 1-D visualizations for each chromosome.

### Running on Slurm

Cactus is a Python script that uses [Toil](https://github.com/DataBiosphere/toil) to execute different programs in parallel. Toil supports distributed computing environments such as [Slurm](https://slurm.schedmd.com/overview.html). Slurm is now the best way to run Cactus at scale. They key challenge is to make sure that each cactus job gets a good memory estimate.  If the memory estimate is too high, then the job will use too much cluster resources or, even worse, fail completely because it asks for more resources than are available on the cluster. If the memory estimate is too low, then Slurm may evict the job for using too much memory.  Both of these errors are very difficult to recover from, unfortunately.  Cactus does its best to estimate the memory from the input data (and should do fine in the current tutorial), but there are three options to override the memory estimates of bigger jobs:

* `--mgMemory` : Memory for minigraph construction.
* `--consMemory` : Memory for cactus alignment.
* `--indexMemory` : Memory for full-genome vg indexing.

See [here](https://github.com/ComparativeGenomicsToolkit/cactus/blob/hprc-v1.1/doc/mc-pangenomes/hprc-v1.1-mc.md#pangenie-filtering) for an example of how Cactus was run on Slurm to generate the v1.1 HPRC graphs. 

To run the previous `cactus-pangenome` command on Slurm instead of locally, you must do the following.

First, install the Cactus virtual environment:

```
wget -q https://github.com/ComparativeGenomicsToolkit/cactus/releases/download/v2.6.13/cactus-bin-v2.6.13.tar.gz
```
Then follow the instructions [described here](https://github.com/ComparativeGenomicsToolkit/cactus/blob/v2.6.13/BIN-INSTALL.md)

Next, switch to the directory where your input data is and where you want to run cactus and run the following. 

**IMPORTANT** Toil/Cactus do not (yet) understand cluser time limits. This will change soon (our cluster will be adopting time limits this month), but in the meantime, you need to make sure that the default time limit for all jobs is longer than the slowest job (which is almost always minigraph construction). One way to do this is with the `TOIL_SLURM_ARGS` environment variable. In general, this variable lets you add any options you want to every job submitted to the cluster by Cactus (see `sbatch --help` for a listing) of possible options. If you do not specify this, jobs will be submitted with some default limit (3 hours, I think) and get evicted if they go longer. Thanks **Mamana Mbiyavanga** for helping to figure this out!!!

**Update:** Previous versions of this document used `--giraffe clip filter` below.  Since Cactus v2.8.5, this is [no longer necessary](https://github.com/ComparativeGenomicsToolkit/cactus/pull/1424): you can use just `--giraffe filter` (or leave the `--giraffe` option out altogether to go all in on haplotype subsampling). Not using `--giraffe clip` drastically reduces peak memory usage. 

```
export TOIL_SLURM_ARGS="-t 1440"

rm -rf slurm-logs ; mkdir -p slurm-logs

cactus-pangenome ./js ./hprc10.seqfile --outDir ./hprc10 --outName hprc10 --reference GRCh38 CHM13 \
--filter 2 --haplo --giraffe filter --viz --odgi --chrom-vg clip filter --chrom-og --gbz clip filter full \
--gfa clip full --vcf --vcfReference GRCh38 CHM13 --logFile ./hprc10.log
--consCores 8 --mgMemory 128Gi --batchSystem slurm --batchLogsDir ./slurm-logs --binariesMode singularity
```

These are the differences with the previous example:
* `--workDir` not set: but you can put a location that's available on all worker nodes if you want.
* `--batchSystem slurm`: activates slurm support
* `--batchLogsDir ./slurm-logs` : improves logging of slurm-specific issues by keeping a local copy of all slurm logs
* `--binariesMode singularity` : run all cactus binaries from inside the singularity container

If you see jobs disappearing or mysteriously stopping and restarting, try running `cat *` in the `--batchLogsDir` directory to see if there are any errors.  For example, if you are losing jobs due to exceeding the time limit (see discussion of `TOIL_SLURM_ARGS` above), then you will see some logs with `CANCELLED DUE TO TIME LIMIT` in them. 

## Part 2: Pangenome Graph Properties

**NOTE**: You may not exactly reproduce the exact numbers below, even running on the same data with the same version as Cactus is not deterministic due to how it is parallelized.  Your numbers should still be extremely close, though.

### Basic Statistics

The very first thing to check is the size of your graph.  You can do this with `vg stats -lz`:

```
singularity exec -H $(pwd) docker://quay.io/comparative-genomics-toolkit/cactus:v2.6.13 \
vg stats -lz ./hprc10/hprc10.gbz

```

It will show the number of nodes and edges in the graph along with the total sequence length over all nodes:
```
nodes	28195250
edges	38322112
length	3145521882
```

You can compare that to the length of GRCh38 in the graph
```
singularity exec -H $(pwd) docker://quay.io/comparative-genomics-toolkit/cactus:v2.6.13 \
vg paths -x ./hprc10/hprc10.gbz -S GRCh38 -E | awk '{sum += $2} END {print sum}'
```
which is `3099922541`. There is `~45Mbp`bp of additional (excluding most heterochromatic) sequence added to the pangenome from CHM13 and the four samples. Something on the order of a few megabases per sample is reasonable.  If your results are much different, then that is a definite warning sign that something went very wrong.

Looking at the `.full` graph shows how much additional sequence is added by the centromeres.

```
singularity exec -H $(pwd) docker://quay.io/comparative-genomics-toolkit/cactus:v2.6.13 \
vg stats -lz ./hprc10/hprc10.full.gbz
```
An extra gigabase in this case. We cannot effectively index or map to such graphs (centromere alignment is something we are actively working on, though!)
```
nodes	29372041
edges	39555335
length	4213877926
```

You can use `vg paths` to inspect the amount of sequence (total length of all embedded paths) of any given sample of haplotype.  For example

```
singularity exec -H $(pwd) docker://quay.io/comparative-genomics-toolkit/cactus:v2.6.13 \
vg paths -x ./hprc10/hprc10.gbz -S HG00438 -E | awk '{sum += $2} END {print sum}'

```
```
singularity exec -H $(pwd) docker://quay.io/comparative-genomics-toolkit/cactus:v2.6.13 \
vg paths -x ./hprc10/hprc10.gbz -Q HG00438#1 -E | awk '{sum += $2} END {print sum}'

```
```
singularity exec -H $(pwd) docker://quay.io/comparative-genomics-toolkit/cactus:v2.6.13 \
vg paths -x ./hprc10/hprc10.gbz -Q HG00438#2 -E | awk '{sum += $2} END {print sum}'
```

Show that there is `5679580423`bp for `HG00438`, with `2841204110` and `2838376313` in its first (paternal) and second (maternal) haplotype, respectively. 

The aforementioned `hprc10/chrom-subproblems/contig_sizes.tsv` gives a breakdown of the length of each haplotype in each chromosome. Can be useful to load into a spreadsheet and/or graph in order to check that all input haplotypes are properly represented in the graph.

`minigraph-cactus` graphs are linearized along the reference genome (GRCh38 in this case).  There is exactly one graph component for each contig in GRCh38.  And each component has exactly two tips or stubs (nodes with zero edges at one of their ends). You can count the tips in the graph with

```
singularity exec -H $(pwd) docker://quay.io/comparative-genomics-toolkit/cactus:v2.6.13 \
vg stats -HT ./hprc10/hprc10.gbz | sed -e 's/heads//g' -e 's/tails//g' | wc -w
```
Giving a result of `390`. This is two times the number of contigs in GRCh38, `195`, which can be inspected with

```
wget -q https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/references/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta.gz
zcat GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta.gz | grep '>' | wc -l
```

To verify the number of graph components, you can use `vg chunk` to break up the graph by chromosome (there isn't really a practical use for this as `cactus-pangenome` will have already output chromosome graphs (see above).

```
mkdir chrom-components
singularity exec -H $(pwd) docker://quay.io/comparative-genomics-toolkit/cactus:v2.6.13 \
vg chunk -x ./hprc10/hprc10.gbz -C -b ./chrom-components/chunk -t 32

```

This will make `195` (one for each GRCh38 contig) `.vg` files in `./chrom-components`

```
ls chrom-components/*.vg | wc -l
195
```

### Subgraph Extraction

Working with whole-genome, or even chromosome, graphs can be unwieldy for many tools such as those used for visualization. You can use `vg` or `odgi` to extract smaller regions from a larger graph.  This tutorial will focus on `vg`, but remember that the `odgi` commands you learned in the PGGB tutorial will apply to `.og` output of `cactus-pangenome` as well.

**Note**: `vg` can read `.gfa`, `.vg`, `.gbz` and `.xg` files, though running times and memory usage can vary substantially.


The simplest way to extract a subgraph is by performing queries on `GRCh38` coordinates using `vg chunk` on the `.gbz` file. For example to extract the `lrc_kir` region for visualization with Bandage-NG (which expects `.gfa`), use

```
singularity exec -H $(pwd) docker://quay.io/comparative-genomics-toolkit/cactus:v2.6.13 \
vg chunk -x ./hprc10.gbz -S ./hprc10.snarls -p GRCh38#0#chr19:54015634-55094318 -O gfa > ./lrc_kir.gfa
```

The command looks for the region `chr19:54015634-55094318` in `GRCh38`, then pulls out the smallest site (aka bubble aka snarl) in the graph that contains it (this is what `-S ./hprc10.snarls` is used for).

Instead of pulling out the site, you can use `-c` to specify the number of steps away from the selected region to extract.  For example, `-c 0` will only extract the given path with no added variation.

Subgraph extraction from `.gbz` is rather slow, and does **not** return non-reference paths in the subgraphs. So in the example abovel `lrc_kir.gfa` will only have paths for CHM13 and GRCh38.  If you want to add in the other haplotypes, you can try adding `-T` to `vg chunk`.

If you will be making many queries and/or you want to query on non-reference genomes, the easiest thing may be to create an `xg` index.  This is how to make an `.xg` index of the full graph.  This one will be more appropriate for querying samples that are not `GRCh38` and therefore will potentially be fragmented in the default graph.  (unfortunately `vg chunk` does not yet transparently handle querying on path fragments).  You can also use the `odgi extract` with the `hprc10.full.og` graph that was already made.

```
singularity exec -H $(pwd) docker://quay.io/comparative-genomics-toolkit/cactus:v2.6.13 \
bash -c "vg convert  ./hprc10.full.gbz -x > ./hprc10.full.xg"
```

(we need to use `bash -c "<command>"` in order for the `>` redirect to work)

### Format Conversion

Use `vg convert` to convert between formats and `vg stats -F` to determine which format a graph is.

The most important formats for `vg` are:

* `.gbz` : Highly compressed haplotype paths (will scale to 1000s of samples). Read-only. Required for `vg giraffe`.
* `.vg` : Less compressed but can be modified (ie with `vg mod`)
* `.xg` : Indexed for fast lookup on any path. Read-only.
* `.gfa` : Standard text-based interchange.  By default, `vg` stores paths a W-lines, but `vg convert` can be used to change to `P`-lines.  `vg` **cannot** read `.gfa.gz`.
* `.og` : ODGI format, which combines best properties of `.xg` and `.vg` but will not scale as well as well `.gbz`. 

In general, you will probably mostly use `.gbz`, 'gfa' and `.og` files for your graphs.

### Panacus

[Panacus](https://github.com/marschall-lab/panacus) is a tool for making beautiful figures describing the coverage in a pangenome graph. It is not (yet) included in the Cactus Docker image, but you can locally install it as follows (see other options in its manual):

With Conda
```
mamba install -c conda-forge -c bioconda panacus
```

or manually with a Python virtualenv:
```
virtualenv -p python3 venv-panacus
. venv-panacus/bin/activate
pip install -U matplotlib numpy pandas scikit-learn scipy seaborn

wget --no-check-certificate -c https://github.com/marschall-lab/panacus/releases/download/0.2.3/panacus-0.2.3_linux_x86_64.tar.gz
tar -xzvf panacus-0.2.3_linux_x86_64.tar.gz
# suggestion: add tool to path in your ~/.bashrc
export PATH="$(readlink -f panacus-0.2.3_linux_x86_64/bin)":$PATH
```

And go through all the examples on the [Panacus](https://github.com/marschall-lab/panacus) webpage using this graph as input.  Note that even examples that use PGGB graphs can still be run on your graph.

**IMPORTANT** To exclude reference paths, use `grep -ive 'grch38\|chm13'` instead of `grep -ve 'grch38\|chm13'` and `grep ^W` instead of `grep ^P`, as well as add `|sort | uniq` when making paths lists. 

So to run the first example from the Panacus website, you would do the following:

```
gzip -d hprc10.gfa.gz

grep '^W' hprc10.gfa | cut -f2 | grep -ive 'grch38\|chm13' | sort | uniq > hprc10.paths.haplotypes.txt

RUST_LOG=info panacus histgrowth -t8 -l 1,2,1,1,1 -q 0,0,1,0.5,0.1 -S -a -s hprc10.paths.haplotypes.txt hprc10.gfa > hprc10.histgrowth.node.tsv

panacus-visualize -e hprc10.histgrowth.node.tsv > hprc10.histgrowth.node.pdf
```

Using `panacus`, especially with `--count bp` to chart sequence length instead of nodes, is a very good way to visualize the diversity of the samples in the graph.

## Part 3: Mapping Reads to the Graph

### Short Read Mapping

The best way to map short read (genomic) data to the pangenome is `vg giraffe`.  `vg giraffe` supports single and paired-end inputs, and will output graph mappings in [GAM](https://github.com/vgteam/vg/wiki/File-Formats#gam-graph-alignment--map-vgs-bam) or [GAF](https://github.com/lh3/gfatools/blob/master/doc/rGFA.md#the-graph-alignment-format-gaf) format.  See the [surjecting](#surjecting-to-bam) section below for details on outputting BAM.

#### Test Reads

You can find some 30X paired-end reads from Genome in a Bottle's HG002 here. Each file is `~33Gb`:
```
wget https://storage.googleapis.com/brain-genomics-public/research/sequencing/fastq/hiseqx/wgs_pcr_free/30x/HG002.hiseqx.pcr-free.30x.R1.fastq.gz
wget https://storage.googleapis.com/brain-genomics-public/research/sequencing/fastq/hiseqx/wgs_pcr_free/30x/HG002.hiseqx.pcr-free.30x.R2.fastq.gz
```

#### Mapping to the Allele-Frequency Filtered Graph (the old way)

You can map the above reads with `giraffe` using this command (it assumes the reads are in the same location as the graph, but you can modify it accordingly, even adding another `-v` argument if necessary):

```
singularity exec -H $(pwd) docker://quay.io/comparative-genomics-toolkit/cactus:v2.6.13 \
bash -c "vg giraffe -Z ./hprc10/hprc10.d2.gbz -f ./hprc10/HG002.hiseqx.pcr-free.30x.R1.fastq.gz -f ./hprc10/HG002.hiseqx.pcr-free.30x.R2.fastq.gz -o gaf | bgzip > ./hprc10/hprc10.hg002.gaf.gz"
```

This takes about 2.5 hours and takes 47 Gb of RAM.

#### Mapping to the Personal Pangenome (the new way -- not yet published)

Instead of mapping to the filtered graph, you can use the reads to extract a personal pangenome and map to that. This way you keep rare variants that are present in the sample, and exclude common ones that aren't. In most benchmarks so far, this helps accuracy of downstream applications.  It does require one extra step, though, which is extracting the kmers from the input reads with `kmc`.

First you need to make a file containing the paths of your reads. Assuming it's called `hg002.reads.txt` in the current directory:

```
printf "./HG002.hiseqx.pcr-free.30x.R1.fastq.gz\n../HG002.hiseqx.pcr-free.30x.R2.fastq.gz\n" > hg002.reads.txt
```

Then you use `kmc` to make the kmers index (`hg002.kff`)

```
singularity exec -H $(pwd) docker://gregorysprenger/kmc:v3.2.2 \
kmc -k29 -m128 -okff -t32 @./hg002.reads.txt hg002 .
```

which takes about 15 minutes and 128Gb of memory.  The `.` as the last argument is telling `kmc` to use the current working directory as its workind directory. 

And you use this index to map to the unfiltered graph with `vg giraffe`. 

```
singularity exec -H $(pwd) docker://quay.io/comparative-genomics-toolkit/cactus:v2.6.13 \
bash -c "vg giraffe -Z ./hprc10/hprc10.gbz -f ./hprc10/HG002.hiseqx.pcr-free.30x.R1.fastq.gz -f ./hprc10/HG002.hiseqx.pcr-free.30x.R2.fastq.gz -o gaf --sample HG002 --progress --kff-name ./hg002.kff --haplotype-name ./hprc10/hprc10.hapl | bgzip > ./hprc10/hprc10.hg002.new.gaf.gz"
```

This takes about 2.75 hours and 64Gb of RAM, and also produces `./hprc10.HG002.gbz`, which is the personal pangenome graph itself.

### Long Read Mapping

Note: this section is adpated from [methods for the mc paper](https://github.com/ComparativeGenomicsToolkit/cactus/tree/master/doc/mc-paper/hprc#long-read-mapping-with-graphaligner).

`vg giraffe` will soon be able to map long reads, but is not ready yet. For now, you should use [GraphAligner](https://github.com/maickrau/GraphAligner).

In order for the resulting mapping to be compatible with the `.gbz`, we first convert the `.gbz` into a `.gfa`.  You can run `GraphAligner` on the `.gfa` that comes out of `cactus-pangenome`, but the [coordinates will be different](https://github.com/ComparativeGenomicsToolkit/cactus/blob/master/doc/pangenome.md#node-chopping) from the `.gbz` (note the `--vg-algorithm` flag is important for ensuring this).

```
singularity exec -H $(pwd) docker://quay.io/comparative-genomics-toolkit/cactus:v2.6.13 \
bash -c "vg convert ./hprc10/hprc10.gbz -f --vg-algorithm > ./hprc10/hprc10.gbz.gfa"
```

This takes 43 seconds and 10Gb RAM.

Now download some public GIAB HiFi reads (or use your own):

```
# log in as anonymous, put anything for password
ftp ftp-trace.ncbi.nlm.nih.gov
cd /giab/ftp/data/AshkenazimTrio/HG002_NA24385_son/PacBio_CCS_15kb_20kb_chemistry2/reads/
get m64011_190830_220126.fastq.gz
```

*Protip*: or use `aria2c` to get the file way faster
```
aria2c -x8 -j8 -s8 ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/HG002_NA24385_son/PacBio_CCS_15kb_20kb_chemistry2/reads/m64011_190830_220126.fastq.gz
```

Now run GraphAligner in pangenome mode. I'm using a relatively old version of `GraphAligner` here since the latest one, `1.0.17b--h21ec9f0_2`, doesn't seem to work properly (tons of assert errors, about 100X slower than expect)

```
singularity exec -H $(pwd) docker://quay.io/biocontainers/graphaligner:1.0.13--he1c1bb9_0  \
GraphAligner -g ./hprc10/hprc10.gbz.gfa -f m64011_190830_220126.fastq.gz -a HG002.hifi.gam -x vg -t 32

```

This takes about 2.5 hours and 220Gb RAM. 

You can now use the `.gam` and `.gbz` together with with `vg pack/call` as described in the [SV Genotyping with vg](#sv-genotyping-with-vg) section (though it should be passed into `vg pack` with `-g`, not `-a` as you would for `.gaf.gz`).  You do not need to further use `hprc10.gbz.gfa`.

**IMPORTANT** : This version of `GraphAligner` doesn't output mapping qualities. So if you use `vg pack` on it, make sure to not use a quality filter, ie use `-Q0`.

If you had mapped to `hprc10.gfa.gz` instead of converting the `.gbz`, you can still use `vg pack/call`, but you would need to run theme on `hprc10.gfa.gz`.  Also note, `vg call` must always be run on the exact same graph as `vg pack`.  Even if the nodes are identical, if you change formats between these two tools, then `vg call` will crash).

### Surjecting To BAM

#### From GAF/GAM to BAM

You can project your read mappings from the graph to a linear reference with `vg surject`.  This will let you output your mappings in BAM format, which can be used with non-pangenome tools like `DeepVariant`, `samtools`, `GATK` etc.

You can project your mappings to any reference path in the graph (as selected with `--reference` in `cactus-pangenome`), so GRCh38 or CHM13 in the example.  You can in theory project reads to any sample in the graph (even non reference samples) but it is a little trickier and not covered here (requires updating the `.gbz` with `vg gbwt`).

You must first create a list of reference paths. Note we're using the `.full` graph here just to make our path lists with the `vg paths` command.  This is important when getting the CHM13 path list, as they will be fragmented otherwise. **Important** you don't want to map or call with the `.full` graph -- it's just for getting this path list. 

```
singularity exec -H $(pwd) docker://quay.io/comparative-genomics-toolkit/cactus:v2.6.13 \
vg paths -x ./hprc10/hprc10.full.gbz -S GRCh38 -L > grch38.paths.txt

singularity exec -H $(pwd) docker://quay.io/comparative-genomics-toolkit/cactus:v2.6.13 \
vg paths -x ./hprc10/hprc10.full.gbz -S CHM13 -L > chm13.paths.txt
```

To project your reads to `GRCh38`, do the following (use `-p chm13.paths.txt` to instead project to CHM13).  If you don't supply a path list with `-F` it will project to to a mix of GRCh38 and CHM13 which is almost certainly *not* what you want.

```
singularity exec -H $(pwd) docker://quay.io/comparative-genomics-toolkit/cactus:v2.6.13 \
bash -c "vg surject -x ./hprc10/hprc10.gbz -G ./hprc10/hprc10.hg002.new.gaf.gz --interleaved -F grch38.paths.txt -b -N HG002 -R 'ID:1 LB:lib1 SM:HG002 PL:illumina PU:unit1' > ./hprc10/hprc10.hg002.new.bam"
```

This takes about 4 hours and 64Gb RAM. 

It's important to use `--interleaved` to tell `surject` that the reads are paired.  The readgroup `-R` is boilerplate tags to help `DeepVariant` or other tools that expect this info in the BAM header.  You can also modify the header yourself with `samtools` if needed.

You should be able to use `hprc10.gbz` for surjection whether you aligned to `hprc10.d2.gbz` or `hprc10.gbz` or the personalized pangenome initially.

#### Mapping Directly to BAM

If you are only ever going to use the BAM, you don't need to create the GAF with `giraffe` then `surject` afterwards -- you can do both at once (use `--ref-paths chm13.paths.txt` to project to CHM13 instead):
```
singularity exec -H $(pwd) docker://quay.io/comparative-genomics-toolkit/cactus:v2.6.13 \
bash -c "vg giraffe -Z ./hprc10/hprc10.gbz -f ./hprc10/HG002.hiseqx.pcr-free.30x.R1.fastq.gz -f ./hprc10/HG002.hiseqx.pcr-free.30x.R2.fastq.gz -o bam --sample HG002 --progress --kff-name ./hg002.kff --haplotype-name ./hprc10/hprc10.hapl -R 'ID:1 LB:lib1 SM:HG002 PL:illumina PU:unit1' --ref-paths ./grch38.paths.txt > ./hprc10/hprc10.hg002.new.bam"

```

This takes about 4 hours and 64 Gb RAM.

See note below about fixing the BAM header if you have surjected onti a different reference than was used for the graph (ie you've surjected onto `CHM13` with a `GRCh38`-based graph.

## Part 4: Genotyping and Variant Calling

### Variant Calling with DeepVariant

[DeepVariant](https://github.com/google/deepvariant) is a state of the art variant caller. It does not use pangenome formats, and rather works on FASTA and BAM files, but has been trained to support data from `vg giraffe / surject`.

First, make a FASTA file from your graph (it is generally best to make the FASTA from the graph, to make sure it matches up exactly.  If you are using a different reference, ie CHM13, use the `.full` graph for this step):

```
singularity exec -H $(pwd) docker://quay.io/comparative-genomics-toolkit/cactus:v2.6.13 \
bash -c "vg paths -x ./hprc10/hprc10.gbz -S GRCh38 -F > ./GRCh38.fa"
singularity exec -H $(pwd) docker://quay.io/comparative-genomics-toolkit/cactus:v2.6.13 \
samtools faidx ./GRCh38.fa
```

**Important** if you are surjecting onto `CHM13` but using the GRCh38-based graph, you need to fix your BAM header as follows:

```
singularity exec -H $(pwd) docker://biocontainers/picard-tools:v2.18.25dfsg-2-deb_cv1 \
PicardCommandLine CreateSequenceDictionary R=./GRCh38.fa O=./GRCh38.dict

singularity exec -H $(pwd) docker://quay.io/comparative-genomics-toolkit/cactus:v2.6.13 \
samtools reheader ./GRCh38.dict ./hprc10/hprc10.hg002.new.bam > ./hprc10/hprc10.hg002.new.fix.bam

mv ./hprc10/hprc10.hg002.new.fix.bam  ./hprc10/hprc10.hg002.new.bam
```

Next, index the BAM (this is a required step for almost any variant caller that takes BAM input)

```
singularity exec -H $(pwd) docker://quay.io/comparative-genomics-toolkit/cactus:v2.6.13 \
samtools sort  ./hprc10/hprc10.hg002.new.bam -O BAM -o ./hprc10/hprc10.hg002.new.sort.bam --threads 8
singularity exec -H $(pwd) docker://quay.io/comparative-genomics-toolkit/cactus:v2.6.13 \
samtools index ./hprc10/hprc10.hg002.new.sort.bam -@ 8
```

These commands took about 30 minutes and very little memory (could be much faster on a local disk).

Finally, run DeepVariant to make a VCF from the BAM.

```
singularity exec -H $(pwd) docker://google/deepvariant:1.6.0 \
  /opt/deepvariant/bin/run_deepvariant \
  --model_type=WGS \
  --ref=./GRCh38.fa \
  --reads=./hprc10/hprc10.hg002.new.sort.bam\
  --output_vcf=./hprc10/hprc10.hg002.new.dv.vcf.gz \
  --output_gvcf=./hprc10/hprc10.hg002.new.dv.g.vcf.gz \
  --make_examples_extra_args="min_mapping_quality=1,keep_legacy_allele_counter_behavior=true,normalize_reads=true" \
  --num_shards=32
```

This took about 13 hours.

### SV Genotyping with vg

We make an important distinction between *genotying* and *calling*:

* *genotyping*: Determine which variants in the graph are present in (each haplotype) of the sample.
* *calling*: Determine which variants in the reads are present in (each haplotype) of the sample. These variants may or may not be in the graph.

One strength of pangenome graphs is that they allow Structural Variants (SVs), which are normally different to determine from short reads, to be efficiently genotyped. One way to do this is with [vg call](https://doi.org/10.1186/s13059-020-1941-7).  This is a two-step process, beginning with a graph alignment (GAF or GAM) from `vg giraffe`.

First, create a `.pack` coverage index (note, if you are using `GraphAligner` output, use `-Q0` below instead):
```
singularity exec -H $(pwd) docker://quay.io/comparative-genomics-toolkit/cactus:v2.6.13 \
vg pack -x ./hprc10/hprc10.gbz -Q5 -a ./hprc10/hprc10.hg002.new.gaf.gz -o ./hprc10/hprc10.hg002.pack
```

This takes about 1 hour and 60 Gb RAM.

Then, create the VCF with `vg call`:

```
singularity exec -H $(pwd) docker://quay.io/comparative-genomics-toolkit/cactus:v2.6.13 \
bash -c "vg call ./hprc10/hprc10.gbz -r ./hprc10/hprc10.snarls -k ./hprc10/hprc10.hg002.pack -s HG002 -S GRCh38 -az | bgzip >  ./hprc10/hprc10.call.vcf.gz"
```

This takes 30 minutes and 40 Gb RAM.

### SV Genotyping with pangenie

**Credit for this section:** [Matteo Ungaro](https://github.com/Overcraft90)

Along with  `vg call`  another tool for Structural Variants detection which works on pangenome graphs is 
[PanGenie](https://github.com/eblerjana/pangenie); however, in contrast with the former  PanGenie  does not leverage any information from reads 
alignment. It uses a HMM to determine the best possible haplotypes combination based on the paths in 
the graph and the short reads dataset for the sample in question.
For this reason,  PanGenie  performs a genome inference for this sample that is represented as a mosaic 
of haplotypes, and which variants are limited those present in the graph space.

#### 1. Installation and Pre-processing

PanGenie  can be installed using some commonly used package manager such as Conda or Singularity. 
Conda installation is the one proposed in this tutorial. First, the tool repository is downloaded with  git 
clone
```
git clone https://github.com/eblerjana/pangenie.git
```
Afterwards, it is possible to access the  `./pangenie`  directory where it is necessary to create a custom 
Conda environment with all the dependencies necessary to run the tool
```
cd pangenie  
conda env create -f environment.yml
```
Finally, the environments has to be activated and tested for the tool’s build
```
conda activate pangenie   
mkdir build; cd build; cmake .. ; make
```

Importantly, the  PanGenie  pipeline later introduced will make use of Snakemake that needs to be present 
as a separate Conda environment, activated when using the tool.

##### 1.2 Filtering and decomposition

The first step to use  PanGenie  is to make sure the graph VCF presents the following properties (more on 
this at the GitHub page for the tool):
* phased samples
* sequence-resolved
* non-overlapping variants
* diploid (important, since in our example CHM13 will be haploid in the VCF)

We need to make sure the VCF is simplified as much as possible and diploid.  To that end,
we start with the VCF output of MC, cut the other reference out (CHM13 in our example),
then normalize it with `vcfwave`. 

```
singularity exec -H $(pwd) docker://ghcr.io/pangenome/pggb:latest \
bash -c "bcftools view -s ^CHM13 ./hprc10/hprc10.vcf.gz \
| vcfwave -I 10000 -t 32 -n | bgzip > hprc10/hprc10.wave.vcf.gz \
&& tabix -fp vcf hprc10/hprc10.wave.vcf.gz"
```

##### 1.3 Running  PanGenie

PanGenie  can be run form command line, but it is recommended to use the pipeline available because the 
output VCF for the sample will be otherwise unviable for conversion to biallelic sites, a post-processing 
step which improves performance. 
The pipeline runs from a  `config.yaml`  file found at this path:
```
/path/to/pangenie/pipelines/run-from-callset/config.yaml)
```

once the tool is installed.

And this is how it looks like:

```
# input vcf with variants to be genotyped (uncompressed)
vcf: /path/to/<decomposed_and_filtered>.vcf
# reference genome
reference: /path/to/<reference>.fna
# path to reads (FASTQ format, uncompressed)
# for each sample that shall be genotyped
reads:
 sample1: /path/to/<sample_name>.fastq
 #sample2: reads-sample2.fastq
# path to PanGenie exectuable
pangenie_genotype: /path/to/pangenie/build/src/PanGenie
pangenie_index: /path/to/pangenie/build/src/PanGenie-index
# name of the output directory (keep in mind PanGenie wants a specific directory where to save all the outputs 
# for instance, <sample_name>-genome_inference
outdir: /path/to/outdir
```

Once you've edited this to put the full path of your reads and `.wave` vcf (**IMPORTANT: Pangenie does not accept gzipped input: you must uncompress everything first**).


Once done,  PanGenie  can be run in interactive from the folder where the  config.yaml  lives or, alternatively, 
form a bash script simply typing
```
snakemake --cores n_of_threads
```
Beware,  PanGenie  requires that the chromosome/contig names match between the reference genome and 
the graph VCF. Normally, it is good practice to keep the reference names unchanged and set the 
chromosome names in the VCF to match the one in the reference; this can be easily done with  sed . For a 
minigraph-CACTUS  pangenome this is how it is done:

```
# example based on a GRCh38 graph
sed -i 's/GRCh38#0#//' <decomposed_and_filtered>.vcf
```

#### 2. Example BASH script

```
#!/bin/bash
#
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=<n_of_threads>
#SBATCH --time=<time_in_hh:mm:ss>
#SBATCH --mem=<allocated_ram>
#
#SBATCH --job-name=pangenie
#SBATCH --output=genome_inference.out
#
#SBATCH --partition=<partition_name>
source /path/to/anaconda3/etc/profile.d/conda.sh
conda activate snakemake
cd /path/to/pangenie/pipelines/run-from-callset/
snakemake --cores n_of_threads
```

PanGenie  is very fast but the process is driven by the graph complexity and the quality/size of the short 
reads dataset for the sample analyzed. A typical run can take ~1.5 hours and up to 80Gb of RAM. 
However, **1.2 Filtering and decomposition** can be a fairly long process mainly because  `vcfwave`  
multithreads only per chromosome and not over the whole VCF file at once.

#### 3. Post-processing

The results are found in a subfolder of the  `/outdir`  called  `/genotypes` . In there, a single VCF file named 
sample1-genotypes.vcf will be present that needs to be converted to biallelic sites.
To do so it is best to use the available script at the GitHub page for the tool (`convert-to-biallelic.py`) which 
can be run as follow:








