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
cactus-align ./jobstore ./examples/evolverPrimates.txt primates.paf primates.hal --pangenome --pafInput --realTimeLogging --outVG --acyclic simChimp
```

## HPRC Chromosome-by-Chromosome Graph Construction

###Input

*Todo: replace with public-facing links once available*

* `seqfile` containing links to gzipped fasta files for 90 haplotypes from `2020AUG26` bucket + hg38 (no alts) + [CHM13](https://s3.amazonaws.com/nanopore-human-wgs/chm13/assemblies/chm13.draft_v1.0.fasta.gz)
* minigraph GFA `GRCh38-freeze1.gfa.gz`

### Startup

Make a node on AWS with `toil launch-cluster --leaderNodeType t2.medium --leaderStorage 32 ...` and set up cactus as [described here](./running-in-aws.md).

### Preprocessing (about 7 hours)

In order to limit spurious alignments and graph complexity in difficult regions, it is recommended to mask these regions out from the get-go.  In an abundance of caution, we use two methods to detect these regions.  First, align the sequences to create a PAF

```
cactus-graphmap <jobstore> seqfile.txt <minigraph GFA> s3://<bucket>/GRCh38-freeze1-orig.paf --realTimeLogging --outputFasta s3://<bucket>/GRCh38-freeze1.gfa.fa --logFile graphmap-1.log --batchSystem mesos --provisioner aws --nodeTypes r3.8xlarge:0.7 --maxNodes 20 --defaultPreemptable --betaInertia 0 --targetTime 1 
```

Next, we combine the coverage gaps in the PAF with dna-brnn to produce the final masking, using a 100kb length threshold (default in cactus config).  Note that we need to create a `seqfile.masked.txt` seqfile specifying the locations of the masked fasta sequences.  

```
cactus-preprocess <jobstore> seqfile.txt seqfile.masked.txt --realTimeLogging --logFile preprocess.log  --batchSystem mesos --provisioner aws --nodeTypes r3.8xlarge:0.7 --maxNodes 25 --defaultPreemptable --betaInertia 0 --targetTime 1  --maskPAF s3://<bucket>/GRCh38-freeze1-orig.paf  --maskAlpha --brnnCores 8
```

Note that instead of softmasking sequences and leaving them unaligned, they can be clipped out entirely using `--clipAlpha` instead of `--maskAlpha` above.

### Map contigs to minigraph (about 2 hours)

Use `cactus-graphmap` to do the mapping.  It will produce a PAF file of pairwise alignments for cactus to build on.  It will be much faster this time than above, as it will ignore the masked sequences.  A future version of Cactus will clean up the interface in order to make this second call unnecessary! 

```
cactus-graphmap <jobstore> `seqfile.masked.txt` <minigraph GFA> s3://<bucket>/GRCh38-freeze1.paf --realTimeLogging --logFile graphmap2.log --batchSystem mesos --provisioner aws --defaultPreemptable --nodeTypes r3.8xlarge:0.7 --maxNodes 20 --outputFasta s3://<bucket>/GRCh38-freeze1.gfa.fa --maskFilter 100000
```

### Split the sequences by reference chromosome (about 2 hours)

Cactus hasn't yet been tested on the entire input set at once.  In the meantime, we split into chromosomes and make a graph for each.  All output will be written in bucket specified by `--outDir`.  A separate cactus project will be written for each chromosome (`--refContigs`).  The decomposition is determined from mapping: an input contig maps to the reference chromosome to which it has the most alignment in the PAF. Input contigs that cannot be assigned confidently to a reference chromosome will be flagged as `_AMBIGUOUS_`. It's currently important to use `--reference hg38` in order for it not to be filtered out with ambiguity filter.

```
cactus-graphmap-split <jobstore> `seqfile.masked.txt` <minigraph GFA> s3://<bucket>/GRCh38-freeze1.paf --refContigs $(for i in $(seq 1 22; echo X; echo Y; echo M); do echo chr$i; done) --reference hg38 --outDir s3://<bucket>/chroms --realTimeLogging --logFile graphmap-split.log --batchSystem mesos --provisioner aws --defaultPreemptable --betaInertia 0 --targetTime 1 --nodeTypes r3.8xlarge:0.7 --maxNodes 20 
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
cactus-align-batch  ./chromfile.txt s3://<bucket>/align-batch --alignCores 32 --alignCoresOverrides chr1,64 chr2,64, chr3,64 chr4,64, chr5,64 --alignOptions "--pafInput --pangenome --outVG --barMaskFilter 100000 --realTimeLogging"  --batchSystem mesos --provisioner aws --defaultPreemptable  --nodeTypes r4.16xlarge:1.7,r4.8xlarge:0.7 --nodeStorage 1000 --maxNodes 5,20 --betaInertia 0 --targetTime 1 --logFile align-batch.log --realTimeLogging
```

### Combining the output into a single graph

The chromosome graphs can then be merged with:
```
vg combine $(for i in $(seq 1 22; echo X; echo Y; echo M); do echo chr${i}.vg; done) > GRCh38-freeze1.cactus.vg
```


