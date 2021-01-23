The Cactus Pangenome Pipeline
===

## Introduction

[Cactus](../README.md) uses a phylogenetic tree as a guide in order to progressively create multiple alignments.  This heuristic allows Cactus scale linearly with the number of input genomes, by decomposiing the work into one alignment per internal (ancestral) node of the tree.  If the guide tree is fully resolved (binary), only two genomes (plus up to three outgroups) are aligned in each subproblem.

Progressively aligning up a guide tree makes sense when the evolution of the input genomes can be explained by a tree.  It is robust to small errors in the tree, as well as small numbers of non-treelike events (ex incomplete lineage sorting or horizontal genen transfer), making it [suitable for alignments of different vertebrate species](https://doi.org/10.1038/s41586-020-2871-y).

But the tree-like assumption breaks down when considering an alignment of individuals from the *same species*.  Such within-population genome alignments are increasingly in demand as high-quality assemblies become more available (ex: [HPP](https://humanpangenome.org/)), given their potential to better identify and represent structural variation than more traditional reference-based re-sequencing approaches.

The Cactus Pangenome Pipeline adapts [Cactus](../README.md) to no longer rely on a guide tree, by taking advantage of the relative similarity of the input sequence to use minimizer sketches to determine initial anchors, then partial order alignments to refine them.  It also provides the options to generate output in standard pangenome graph formats such as [vg](https://github.com/vgteam/vg) and [GFA](https://github.com/GFA-spec/GFA-spec), in addition to the usual HAL. 

*This is a work in progress and is not yet published.* 

## Overview

The interface is very similar to the [default utilization of Cactus](../README.md), and is dependent on a [seqfile mapping genome names to fasta locations](seqFile-the-input-file).  The main difference here is that a tree need not be provided at the top of the seqfile.  If a tree is present, it must be a star tree (all leaves connected to one root).

Also, a [minigraph](https://github.com/lh3/minigraph) GFA is required.  It can be constructed by running `minigraph -xggs`.  It is suggested but not required to use the fasta files from the seqfile as input here.  Note that the first sequence passed to minigraph will be considered the "reference", and its paths will be acyclic in the graph output.  

The following two Cactus commands are run to produce an alignment and pangenome graph from the minigraph GFA and fasta input:

1. `cactus-graphmap`: Align each input fasta sequence to the minigraph
2. `cactus-align --pangenome --pafInput --outVG`: Run cactus in pangenome mode to produce a HAL alignment and vg graph from the minigraph alignments.

The pipeline can be run without the minigraph GFA by using `cactus-refmap` (which will require identifying one of the inputs as the reference) instead of `cactus-graphmap`, but this will at some cost to sensitivity.

While lastz preprocessing with `cactus-preprocess` is strongly recommended for the Progressive Cactus pipeline, it is not necessary here.  But masking centromeres with the new `cactus-preprocess --alphaMask` can be useful in some cases.  

`cactus-align` does not yet scale to whole human genomes, but this can be worked around by decomposing into chromosomes as [described in the example below](hprc-graph-construction).

## Evolver Simulation Example

This is a very small example along the same lines as Cactus's "evolverMammals" test.  Begin by constructing the minigraph:
```
# download the fasta sequences
for i in `awk '{print $2}' examples/evolverPrimates.txt; do wget $i; done
# make the minigraph with human as the reference
minigraph -xggs simHuman.chr6 simChimp.chr6 simGorilla.chr6 simOrang.chr6 > primates.gfa
```

Align the sequences back to the graph.  Note that the `--outFasta` option is required.  It will be used to update the seqfile with an entry for the minigraph node sequences.  
```
cactus-graphmap ./jobstore ./examples/evolverPrimates.txt primates.gfa primates.paf --outputFasta primates.gfa.fa --realTimeLogging
```

Create the cactus alignment from the seqfile and PAF, and export the output to vg.
```
cactus-align ./jobstore ./examples/evolverPrimates.txt primates.paf primates.hal --pangenome --pafInput --realTimeLogging --outVG
```

## HPRC Chromosome-by-Chromosome Graph Construction

###Input

*Todo: replace with public-facing links once available*

* `seqfile` containing links to gzipped fasta files for 90 haplotypes from `2020AUG26` bucket + hg38 (no alts) + [CHM13](https://s3.amazonaws.com/nanopore-human-wgs/chm13/assemblies/chm13.draft_v1.0.fasta.gz)
* minigraph GFA `GRCh38-freeze1.gfa.gz`

### Startup

Make a node on AWS with `toil launch-cluster --leaderNodeType t2.medium --leaderStorage 32 ...` and set up cactus as [described here](./running-in-aws.md).

### Preprocessing (about 7 hours)

minigraph is unable to map within centromeres.  While `cactus-bar` can fill in alignments within these regions to some extent, we currently recommend removing them instead.  This substantially reduces complexity in the resulting pangenome, as well as lowers the probability that assembly or alignment errors get passed off as real variation.  `cactus-preprocess` can be used to remove centromeres as follows. 

Make a `config.xml` file where `cpu` attribute of the `dna-brnn` `preprocessorJob` is increased from `2` to `16` (it takes a little under 2 hours to mask a genome with 16 cores).  Also change the `active` attribute of the `checkUniqueHeaders` job from `1` to `0` (to save time). (leaving everything else default is fine)

Make a seqfile for the output fastas in `./masked`
```
cactus-prepare ./seqfile.txt --outDir ./masked > /dev/null
# edit masked/seqfile.txt to remove Anc0, the tree, and verify paths (and put them in S3!!)
```

Run the preprocessor, clipping out brnn-masked regions > 100kb (except on reference). (bed files of what was clipped also written). 
```
cactus-preprocess <jobstore> ./seqfile.txt ./masked/seqfile.txt --clipAlpha 100000 --configFile ./config.xml --ignore hg38 --realTimeLogging --logFile preprocess.log --batchSystem mesos --provisioner aws --defaultPreemptable --nodeTypes r3.8xlarge:0.7 --maxNodes 20
```

Note: To mask the centromeres and leave them as unaligned bubbles in the graph (instead of clipping them out), do the following:
* Pass `--maskAlpha 100000` instead of `--clipAlpha 100000` to `cactus-preprocess` above
* Pass `--maskFilter 100000` to `cactus-graphmap` below to tell `minigraph` to ignore these regions
* Pass `--barMaskFilter 100000` to `cactus-align` below to tell `cactus-bar` to ignore these regions

### Map contigs to minigraph (about 3 hours)

Use `cactus-graphmap` to do the mapping.  It will produce a PAF file of pairwise alignments for cactus to build on.  It will also edit the input seqfile in order to add the minigraph fasta (specfied with `--outputFasta`).  Note that `--refFromGFA hg38` is used in order to not map hg38 and instead pull its alignments from the rGFA tags.  

```
cactus-graphmap <jobstore> `masked/seqfile.txt` <minigraph GFA> s3://<bucket>/GRCh38-freeze1.paf --realTimeLogging --logFile graphmap.log --batchSystem mesos --provisioner aws --defaultPreemptable --nodeTypes r3.8xlarge:0.7 --maxNodes 20 --outputFasta s3://<bucket>/GRCh38-freeze1.gfa.fa --refFromGFA hg38
```

### Split the sequences by reference chromosome (about 1 hour)

Cactus hasn't yet been tested on the entire input set at once.  In the meantime, we split into chromosomes and make a graph for each.  The GFA results can be simply catted to make a single pangenome graph.  All output will be written in bucket specified by `--outDir`.  A separate cactus project will be written for each chromosome (`--refContigs`).  The decomposition is determined from mapping: an input contig maps to the reference chromosome to which it has the most alignment in the PAF. Input contigs that cannot be assigned confidently to a reference chromosome will be flagged as `_AMBIGUOUS_`. It's currently important to use `--reference hg38` in order for it not to be filtered out with ambiguity filter.

```
cactus-graphmap-split <jobstore> `masked/seqfile.txt` <minigraph GFA> s3://<bucket>/GRCh38-freeze1.paf --refContigs $(for i in $(seq 1 22; echo X; echo Y; echo M); do echo chr$i; done) --reference hg38 --outDir s3://<bucket>/chroms --realTimeLogging --logFile graphmap-split.log --batchSystem mesos --provisioner aws --defaultPreemptable --nodeTypes r3.8xlarge:0.7 --maxNodes 20 
```

### Make a pangenome graph for each chromosome (about 20 hours)

Begin by downloading the chromfile and seqfiles that were created above (they must be input from local folders)
```
pip install awscli
aws s3 sync s3://<bucket>/chroms/seqfiles ./seqfiles
aws s3 cp s3://<bucket>/chroms/chromfile.txt ./
```

Use `cactus-align-batch` to align each chromosome on its own aws instance.  The `--alignCoresOverrides` option is used to ensure that the larger chromosomes get run on bigger instances.

```
cactus-align-batch  ./chromfile.txt s3://<bucket>/align-batch --alignCores 32 --alignCoresOverrides chr1,64 chr2,64, chr3,64 chr4,64, chr5,64 --alignOptions "--pafInput --pangenome --outVG --realTimeLogging"  --batchSystem mesos --provisioner aws --defaultPreemptable  --nodeTypes r4.16xlarge:1.7,r4.8xlarge:0.7 --nodeStorage 1000 --maxNodes 5,20 --betaInertia 0 --targetTime 1 --logFile align-batch.log --realTimeLogging
```

### Combining the output into a single graph

It is important to first make sure that the chromosome graphs have unique ids:
```
vg ids -j $(for i in $(seq 1 22; echo X; echo Y; echo M); do echo chr${i}.vg; done)
```

The can then be merged with:
```
vg combine $(for i in $(seq 1 22; echo X; echo Y; echo M); do echo chr${i}.vg; done) > GRCh38-freeze1.cactus.vg
```
