# Minigraph-Cactus Pangenome Construction and Downstream Applications

## Table of Contents

* [Abstract](#abstract)
* [Key Reference Material](#key-reference-material)
* [Part 1: Pangenome Graph Construction](#part-1-pangenome-graph-construction)
     * [Cactus Setup][#cactus-setup]
     * [Input Data](#input-data]
     * [Build and Index the Pangenome Graph](#build-and-index-the-pangenome-graph)

## Abstract

This is a tutorial written to support the **Reference Graph Pangenome Data Analysis Hackathon 2023 Nov. 13-17 in Cape Town, South Africa**. The aim is to provide detailed instructions on how to create a pangenome reference graph with Minigraph-Cactus then use it for some downstream analysis like variant calling and genotyping.

Unlike some previous workshops, and most of the existing Cactus documentation, this tutorial will focus on whole-genome human data. As such, it will need to be run over a period of time longer than a typical workshop session. The running times and memory usage of each command will be given. 

## Key Reference Material

Please visit these links for related material and background information before proceeding further. **The first link is essential and should absolutely be consulted before continuing and the rest are highly recommended.** 

* [Minigraph-Cactus Manual](https://github.com/ComparativeGenomicsToolkit/cactus/blob/master/doc/pangenome.md): This is essential to read, and includes several small examples (with data) that should be run before tackling whole-genomes.
* [Minigraph-Cactus Paper](https://doi.org/10.1038/s41587-023-01793-w): The methods are described in detail here.
* [HPRC v1.1 Minigraph-Cactus Instructions](https://github.com/ComparativeGenomicsToolkit/cactus/blob/master/doc/mc-pangenomes/hprc-v1.1-mc.md): Exact commands and explanations in order to exactly reproduce the latest released HPRC graph. The commands themselves assume a SLURM cluster but can be trivially modified to run on a single computer (remove `--batchSystem slurm`). 
* [HPRC Paper](https://doi.org/10.1038/s41586-023-05896-x): Detailed analysis of the HPRC graph, and examples of many downstream applications of Minigraph-Cactus pangenomes. 
* @jeizennga's [2023 Memphis Workshop](https://github.com/pangenome/MemPanG23/blob/main/lessons/Day_3a_vg_mapping_and_calling.md), which served as an inspiration for this tutorial.

## Part 1: Pangenome Graph Construction

### Cactus Setup

**Important:** We will be using [Cactus v2.6.9](https://github.com/ComparativeGenomicsToolkit/cactus/releases/tag/v2.6.9) for this tutorial. Be warned that it may not work for newer or older versions.

For simplicity, all cactus will be run in "single-machine" mode via its [docker](https://www.docker.com/) image.  Cactus supports distributed computing environments via [slurm](https://github.com/ComparativeGenomicsToolkit/cactus/blob/master/doc/progressive.md#running-on-a-cluster) and [AWS/Mesos](https://github.com/ComparativeGenomicsToolkit/cactus/blob/master/doc/running-in-aws.md).

In order to make sure docker is working, try running the following and verify that you do not get an error. 
```
docker run hello-world
```

You can then pull the Cactus image onto your coputer
```
docker pull quay.io/comparative-genomics-toolkit/cactus:v2.6.9
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

**If you are making a pangenome graph with your own data, this should be the only part you need to change, but do see the sxplanation of the options below as some may require adjustments for different data sizes**. Also, nothing changes if you want to use haploid assemblies -- just do not use the `.1` and `.2` suffixes (see `CHM13` and `GRCh38` above).

### Build and Index the Pangenome Graph

I am going to run on 32-cores in order to simulate my understanding of an "average" node on your cluster. As you've seen in the [Minigraph-Cactus Manual](https://github.com/ComparativeGenomicsToolkit/cactus/blob/master/doc/pangenome.md) (please go back and read it if you haven't already), the simplest way to build the graph is with the `cactus-pangenome` command.

Here it is, with an explanation of each option following below.

```
docker run -it --rm -v $(pwd):/data --user $UID:$GID quay.io/comparative-genomics-toolkit/cactus:v2.6.9 \
cactus-pangenome /data/js /data/hprc10.seqfile --outDir /data/hprc10 --outName hprc10 --reference GRCh38 CHM13 \
--filter 2 --haplo --giraffe clip filter --viz --odgi --chrom-vg clip filter --chrom-og --gbz clip filter full \
--gfa clip full --vcf --vcfReference GRCh38 CHM13 --logFile /data/hprc10.log \
--consCores 8
```

For `docker run`:
* `-it`: interactive / tty.  Boilerplate for most command line tools
* `--rm`: save space by removing the container after running
* `-v $(pwd):/data`: mount the current directory to `/data` in the container
* `--user $UID:$GID`: run as the current user in the container
* `quay.io/comparative-genomics-toolkit/cactus:v2.6.9`: the cactus docker image

For `cactus-pangenome`:
* `/data/js`: Scratch directory that will be created for Toil's jobstore
* `/data/hprc10.seqfile`: The input samples and assemblies. This file was created above.
* `--outDir /data/hprc10`: The output directory. All results will be here. Remember anything relative to `/data` in the docker container will end up in your current working directory that you're running `docker run` from.
* `--outName hprc10`: This will be the prefix of all output files.
* `--reference GRCh38 CHM13`: Specify these two samples as reference genomes. Reference samples are indexed a little different in `vg` to make their coordinates easier to use. Also, the first reference given (GRCh38 in this case), is used to anchor the entire graph and is treated differently than the other samples.  Please see [here for more details](https://github.com/ComparativeGenomicsToolkit/cactus/blob/master/doc/pangenome.md#reference-sample).
* `--filter 2`: Create an Allele-Frequency filtered (AF) graph that contains only nodes and edges supported by at least 2 haplotypes. This can lead to better mapping performance. `--filter 9` was used for the 90-assembly HPRC graph.  
* `--haplo`: We are actually phasing out the Allele-Frequency filtering as described above in favour or dynamic creation of personal pangenomes. Using this option will create the necessary indexes for this functionality.
* `--giraffe clip filter`: Make giraffe indexes for the Allele-Frequency filtered graph in addition to the (default) clipped graph.
* `--viz`: Make an ODGI 1D visualization image for each chromosome.
* `--odgi`: Make an ODGI formatted whole-genome graph
* `--chrom-vg clip filter`: Make VG formatted chromosome graphs for the both the AF filtered and (default) clipped pangenome.
* `--gbz clip filter full`: Make GBZ formatted whole-genome graphs for AF filtered, (default) clipped and full (containing unaligned centromeres) graphs.
* `--gfa clip full`: Make GFA formatted whole-genome graphs for (default) clipped and full graphs.
* `--vcf`: Make a VCF (based on the first reference) version of the graph
* `--vcfReference GRCh38 CHM13`: Specify that we want two VCFs, one for each reference
* `--logFile /data/hprc10.log`: All logging information will end up here in addition to `stderr`.  Important to save!
* `--consCores 4`: Specify 4 threads for each core cactus job (`cactus_consolidated`). By default it will use all cores available on your system.  By reducing to `8`, we attempt to run up to 4 chromosomes at once (assuming 32 cores total).  Note that this will increase peak memory usage.

All of the above is explained in more detail in the [Minigraph-Cactus Manual](https://github.com/ComparativeGenomicsToolkit/cactus/blob/master/doc/pangenome.md). We are erring on the side of producing lots of different indexes, but it's usually easier than going back and regenerating any forgotten ones.

Here are some details about the resources used. I'm running on a shared server on a slow network drive, so your times may be faster.

1) Minigraph Construction :
2) Minigraph Mapping :
3) Cactus Alignment :
 * chr1
