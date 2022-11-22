# The Minigraph-Cactus Pangenome Pipeline

Minigraph-Cactus is included in the [Cactus Software Package](../README.md) and is suitable for aligning similar samples, such as those from the *same species*. See [Progressive Cactus](./progressive.md) for aligning different species. 

Please cite the [Minigraph-Cactus paper](https://doi.org/10.1101/2022.10.06.511217 ) when using Minigraph-Cactus.

## Table of Contents

* [Pangenome Data](#pangenome-data)
* [Quick-Start](#quick-start)
* [Introduction](#introduction)
* [Interface](#interface)
* [Output](#output)
* [Yeast Graph](#yeast-graph)
* [Human Graph](#hprc-graph)
* [Frequently Asked Questions](#frequently-asked-questions)

## Pangenome Data

Some pangenomes constructed with Minigraph-Cactus, along with all material to reproduce, can be found [here](./mc-pangenomes/README.md). 

## Quick-Start

Set up an output directory, and copy over the example seqfile (as it will be modified)
```
mkdir -p primates-pg
cp examples/evolverPrimates.txt primates-pg/evolverPrimates.pg.txt
```

Make the SV graph with `minigraph`
```
cactus-minigraph ./jobstore primates-pg/evolverPrimates.pg.txt primates-pg/primates.sv.gfa.gz  --reference simChimp
```

Make the assembly-to-graph alignments with `minigraph`
```
cactus-graphmap ./jobstore primates-pg/evolverPrimates.pg.txt primates-pg/primates.sv.gfa.gz primates-pg/primates.paf  --reference simChimp --outputFasta primates-pg/primates.sv.gfa.fa.gz
```

Create the Cactus base alignment and "raw" pangenome graph
```
cactus-align ./jobstore primates-pg/evolverPrimates.pg.txt primates-pg/primates.paf primates-pg/primates.hal --pangenome --outVG --reference simChimp 
```

Create and index the final pangenome graph and produce VCF files and vg-giraffe indexes. 
```
cactus-graphmap-join ./jobstore --vg primates-pg/primates.vg --outDir ./primates-pg --outName primates-pg --reference simChimp --vcf --giraffe
```

If it worked properly, the input sequences should show up as contigs in the GFA:
```
zcat primates-pg/primates-pg.gfa.gz | grep '^W' | awk '{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 }' | more
W	simChimp	0	simChimp.chr6	0	596350
W	simGorilla	0	simGorilla.chr6	0	599081
W	simHuman	0	simHuman.chr6	0	597871
W	simOrang	0	simOrang.chr6	0	591073
```

The reference path, which will be treated differently by some vg tools (for efficiency), will be identified in the header:
```
zcat primates-pg/primates-pg.gfa.gz | head -1
H	VN:Z:1.1	RS:Z:simChimp
```

The VCF will be based on the reference path (simChimp) and have a sample for each haplotype :
```
gzip -dc primates-pg/primates-pg.vcf.gz | grep CHROM -A 1
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	simGorilla	simHuman	simOrang
simChimp.chr6	15	>1>4	T	G	60.0	.	AC=1;AF=0.333333;AN=3;AT=>1>2>4,>1>3>4;NS=3;LV=0	GT	0	0	1
```

The three input files needed for giraffe are produced:
```
ls -hs primates-pg/primates-pg.d2.*
4.8M primates-pg/primates-pg.d2.dist  2.5M primates-pg/primates-pg.d2.gbz  6.1M primates-pg/primates-pg.d2.min
```

By default, `--giraffe` will produce frequency filtered indexes, with a default minimum coverage of 2 (hence the `.d2`). This means only nodes covered by two haplotypes will appear in the index. This helps `vg giraffe` performance considerably (though a version of Giraffe that no longer needs it is under development). The dataset here is too small for this to be useful.  To index the clipped but unfiltered graph, use `--giraffe clip` or use `--giraffe full` to index the full, unclipped graph. See more detailed explanations below. 


## Introduction

Minigraph-Cactus uses [minigraph](https://github.com/lh3/minigraph) to construct a pangenome graph of structural variation in a set of input assemblies. The assemblies are then mapped back to this graph using minigraph.  These mappings are used as input to [Cactus](../README.md) to construct a new graph that contains variants of all sizes, allowing the input assemblies to be encoded as embedded paths in the graph. The graph is output in  pangenome graph formats such as [vg](https://github.com/vgteam/vg) and [GFA](https://github.com/GFA-spec/GFA-spec), in addition to the usual [HAL](https://github.com/ComparativeGenomicsToolkit/hal).

Pangenomes from Minigraph-Cactus are indexable for and ready for read mapping with [vg Giraffe](https://github.com/vgteam/vg/wiki/Mapping-short-reads-with-Giraffe). 

Unlike Progressive Cactus, Minigraph-Cactus does depend on a predetermined reference genome.  This genome is guaranteed to be acyclic and unclipped in the final graph.  Other genomes in the graph can still be used as reference coordinate systems, however.  For example, we achieved great variant calling performance when projecting mapped reads on a CHM13-referenced graph onto GRCh38. 

## Interface

### Sample Names

The input is a two-column `seqFile` mapping sample names to fasta paths (gzipped fastas are supported). The seqfile is the same as Progressive Cactus, except a tree is not specified.  Cactus may add a star-tree to this file, but it can be ignored.



**A naming convention must be followed for sample names**: The "." character is used to specify haplotype, and should be avoided in sample names unless it is being used that way.  For haploid samples, just don't use a "`.`".  For diploid or polyploid samples, use the form `SAMPLE.HAPLOTYPE` where `HAPLOTYPE` is `1` or `2` for a diploid sample etc:

```
# Diploid sample:
HG002.1  ./HG002.paternal.fa
HG002.2  ./HG002.maternal.fa

# Haploid sample:
CHM13  ./chm13.fa
```

### Pipeline

1) `cactus-minigraph <jobStore> <seqFile> <outputGFA> --reference`: Construct a minigraph in GFA format (may be gzipped) from a set of FASTA files (may also be gzipped).  This is a very thin wrapper over `minigraph -cxggs`.  The reference is added first and the remainder of samples are added in decreasing order of size.  Use the `--mapCores` option to specify the number of cores.

2) `cactus-graphmap <jobStore> <seqFile> <inputGFA> <outputPAF> --reference`: Map each input assembly back to the graph using `minigraph`.  The number of cores for each mapping job can be set with `--mapCores`.  

3) **(Optional)** `cactus-graphmap-split <jobStore> <seqFile> <inputGFA> <inputPAF> --reference --outDir`: Split the input assemblies and PAF into chromosomes using the rGFA tags in the GFA. Doing so reduces the memory requirements in the following steps.  It assigns each contig to a single chromosome according to the alignment in the input PAF, so all inter-chromosomal events will be filtered out.  Contigs that can't be assigned to a chromosome are deemed "ambiguous" and not considered in later steps.

4) `cactus-align <jobStore> <seqFile> <inputPAF> <outHal> --reference --pangenome --outVG --maxLen`: Compute the Cactus multiple genome alignment from the assembly-to-graph minigraph mappings. The `--maxLen` parameter specifies the maximum gap between minigraph mappings that Cactus will attempt to fill at once, and is recommended to be set to 10000.  If `cactus-graphmap-split` was used, the `cactus-align-batch` interface should be used instead (see examples below).

5) `cactus-graphmap-join <jobStore>  --vg --outDir --outName --reference`: Produce the final graph and indexes. This should be run whether or not `cactus-graphmap-split` was used.

### Clipping, Filtering and Indexing

`cactus-graphmap-join` merges chromosome graphs created by `cactus-align-batch`, and also normalizes, clips and filters the graph in addition to producing some useful indexes.  It can produce up to three graphs (now in a single invocation), and a variety of indexes for any combination of them. The three graphs are the

* `full` graph: This graph is normalized, but no sequence is removed. It and its indexes will have `.full` in their filenames. 
* `clip` graph: This is the default graph. Stretches of sequence `>10kb` that were not aligned to the underlying SV/minigraph are removed.
* `filter` graph: This graph is made by removing nodes covered by fewer than 2 haplotypes from the `clip` graph.  It and its indexes will have `.d2` in their filenames.

The `clip` graph is a subgraph of the `full` graph and the `filter` graph is a subgraph of the `clip` graph. Put another way, any node in the `filter` graph exists with the exact same ID and sequence in the `clip` graph, etc. 

The different graphs have different uses. For instance, the current version of `vg giraffe` performs best on the filtered graph (this will hopefully be soon remedied in an update to vg).  For the HPRC v1.0 graph and paper, we used `d9`. When you pass `--giraffe` to `cactus-graphmap-join`, it will make the giraffe indexes on the filtered graph by default.  But you can override this behaviour to produces the indexes for any of the graphs by passing in any combination of [`full`, `clip` and `filter`] to the `--giraffe` options. For example:

`--giraffe`: Make the giraffe indexes for the filtered graph (default choice).

`--giraffe clip`: Make the giraffe indexes for the clipped graph.

`--giraffe clip filter`: Make the giraffe indexes for both the clipped and filtered graph.

The same type of interface applies to all the output specification options: `--vcf`, `--gbz`, `--gfa`, `--giraffe`, `--chrom-vg`. They can all be used without arguments to apply to the default graph (generally the `clip` graph for everything except `--giraffe` which defaults to the `filter` graph), or with any combination of `full`, `clip` and `filter` to be applied to different graphs.

Note that by default, only GFA is output, so the above options need to be used to toggle on any other output types. 

Different clipping and filtering thresholds can be specified using the `--clip` and `--filter` options, respectively. For larger graphs, you probably want to use `--filter N` where `N` represents about 10% of the haplotypes.  It is indeed a shame to remove rarer variants before mapping, but is a necessity to get the best performance out of (the current version) of `vg giraffe`.  

The `--vcf` option will produce two VCFs for each selected graph type. One VCF is a "raw" VCF which contains nested variants, indicated by the `LV` and `PS` tags. The second VCF is one that has gone through [vcfbub](https://github.com/pangenome/vcfbub) to remove nested sites, as well as those greater than 100kb.  Unless you want to explicitly handle nested variants, you are probably best to use the `vcfbub` VCF.  Switch off `vcfbub` with `--vcfbub 0` or specify a different threshold with `--vcfbub N`.

If you want to use the HAL output, `cactus-graphmap-join` can also merge HAL chromosomes from `cactus-align-batch` with the `--hal` option.  These will never be filtered or otherwise processed.

When merging hal files with `--hal`, it is best to set `--indexCores` such that one core is free for hal merging.  So on a 32-core system, use `--indexCores 31`.  This way the slower indexing jobs can be done in parallel with the also slow hal merging (which itself is single-core).

### Output

* `hal`: Cactus's [native alignment format](./progressive.md#using-the-hal-output) can be used to convert to MAF, build assembly hubs, run liftover and comparative annotation.
* `gfa`: A [standard text-based graph format](https://github.com/GFA-spec/GFA-spec/blob/master/GFA1.md). Minigraph-Cactus uses GFA 1.1 as it represents haplotypes as [Walks](https://github.com/GFA-spec/GFA-spec/blob/master/GFA1.md#w-walk-line-since-v11). You can use `vg convert -gfW` to convert from GFA 1.1 to 1.0 and `vg convert -gf` to convert from 1.0 to 1.1.
* `vcf`: A [standard text-based format](https://en.wikipedia.org/wiki/Variant_Call_Format) that represents a pangenome graph as sites of variation along a reference. VCFs exported from the graph are nested, and by default `vcfbub` is used to flatten them.
* `vg`: [vg](https://github.com/vgteam/vg)'s native packed-graph format, can be read and written by vg but does not scale well with the number of paths.
* `gbz`: A read-only [format that scales extremely efficiently with the number of paths](https://github.com/jltsiren/gbwtgraph/blob/master/SERIALIZATION.md). Readable by `vg` tools and required for `giraffe`.
* `snarls`: The start and end nodes of the bubbles in the graph, as well as their nesting relationships.  Used by some `vg` tools like `call` and `deconstruct`.
* `dist`: Snarl distance index required for `vg giraffe`.
* `min`: Minimizer index required for `vg giraffe`.
* `stats.tgz`: Some stats about how much sequence was clipped, including a BED file of the removed sequence.

## Yeast Graph

This is a small test case whose input data is included in cactus that illustrates how to split by chromosome. 

### Yeast: Getting Started

Below is an example of creating a yeast pangenome chromosome by chromosome, referenced on S288C.  

```
# make the seqfile
mkdir -p yeast-pg
cp ./examples/yeastPangenome.txt yeast-pg/

# make the minigraph
cactus-minigraph ./jobstore  ./yeast-pg/yeastPangenome.txt ./yeast-pg/yeast.sv.gfa  --reference S288C

# map back to the minigraph
cactus-graphmap ./jobstore ./yeast-pg/yeastPangenome.txt ./yeast-pg/yeast.sv.gfa ./yeast-pg/yeast.paf --outputFasta ./yeast-pg/yeast.sv.gfa.fa  --reference S288C
```

### Yeast: Splitting By Chromosome
Now the PAF and GFA `minigraph` output can be used to partition the graph and mappings based on the reference genome's (S288C's) chromosomes:

```
cactus-graphmap-split ./jobstore ./yeast-pg/yeastPangenome.txt ./yeast-pg/yeast.sv.gfa ./yeast-pg/yeast.paf --outDir yeast-pg/chroms  --reference S288C
```

This command makes a cactus subproblem for each reference chromosome.  By default, it uses all contigs in the reference.  A subset can be specified using the `--refContigs` option.

In this example, for instance, the `chrI` data can be found as follows.  This is everything required to run `cactus-align` on it as described previously.

```
ls -hs yeast-pg/chroms/chrI/* yeast-pg/chroms/seqfiles/chrI.seqfile 
264K yeast-pg/chroms/chrI/chrI.gfa   44K yeast-pg/chroms/chrI/chrI.paf  4.0K yeast-pg/chroms/seqfiles/chrI.seqfile

yeast-pg/chroms/chrI/fasta:
total 656K
 68K DBVPG6044.0_chrI.fa.gz   68K S288C_chrI.fa.gz   68K UWOPS034614.0_chrI.fa.gz   72K YPS128.0_chrI.fa.gz
244K _MINIGRAPH__chrI.fa      72K SK1.0_chrI.fa.gz   64K Y12.0_chrI.fa.gz
```

Some contigs cannot be assigned to a reference chromosome.  These end up in the `_AMBIGUOUS_` directory:
```
ls -hs yeast-pg/chroms/_AMBIGUOUS_/*
188K yeast-pg/chroms/_AMBIGUOUS_/_AMBIGUOUS_.paf

yeast-pg/chroms/_AMBIGUOUS_/fasta:
total 1.2M
   0 DBVPG6044.0__AMBIGUOUS_.fa.gz     0 S288C__AMBIGUOUS_.fa.gz  1.2M UWOPS034614.0__AMBIGUOUS_.fa.gz     0 YPS128.0__AMBIGUOUS_.fa.gz
   0 _MINIGRAPH___AMBIGUOUS_.fa        0 SK1.0__AMBIGUOUS_.fa.gz     0 Y12.0__AMBIGUOUS_.fa.gz
```

Here we can see that a few contigs from UWOPS034614 were left unplaced (and would be left out of any future cactus jobs).

```
zcat yeast-pg/chroms/_AMBIGUOUS_/fasta/UWOPS034614__AMBIGUOUS_.fa.gz | grep '>'
>id=UWOPS034614|chrXI
>id=UWOPS034614|chrX
>id=UWOPS034614|chrVII
>id=UWOPS034614|chrVIII
```

The reason why these contigs are unassigned to a chromosome can normally be found in `minigraph.split.log`:

```
Query contig is ambiguous: id=UWOPS034614|chrXI  len=792116 cov=0.573045 (vs 0.5) uf=1.47861 (vs 2)
 Reference contig mappings:
  chrVII: 306989
  chrXI: 453918
--
Query contig is ambiguous: id=UWOPS034614|chrVIII  len=738767 cov=0.481758 (vs 0.5) uf=1.06071 (vs 2)
 Reference contig mappings:
  chrVII: 355907
  chrVIII: 335536
--
Query contig is ambiguous: id=UWOPS034614|chrVII  len=632616 cov=0.407576 (vs 0.5) uf= infinity (vs 2)
Assigned contig to chrXI: id=DBVPG6044|chrXI  len=695907 cov=0.972054 (vs 0.5) uf= infinity (vs 2)
Query contig is ambiguous: id=UWOPS034614|chrX  len=1092164 cov=0.49082 (vs 0.25) uf=1.04957 (vs 2)
 Reference contig mappings:
  chrX: 510741
  chrXIII: 536056

```

It is saying that not enough bases in these contigs aligned to a single reference chromosome, given the 50% threshold and 2X uniqueness factor.  These thresholds are explained, and can be adjusted in, the `<graphmap-split>` section of the [cactus config](../src/cactus/cactus_progressive_config.xml).

The sequence-to-graph PAF files can be refined by adding the `--base` option to `cactus-graphmap`.  This can often improve the accuracy of `cactus-graphmap-split`.  This option will be used in the HPRC example below. 

### Yeast: Batch Aligning the Chromosomes

`cactus-align` can be run individually on each chromosome using the seqFiles created above.  This can be automated using the `cactus-align-batch` script, which takes as input a "chromFile", which is just a list of seqFiles and PAFs.  Such a chromFile was generated by `cactus-graphmap-split` above.

The options we would normally pass directly to `cactus-align` must be quoted and passed via `--alignOptions` here:

```
cactus-align-batch ./jobstore ./yeast-pg/chroms/chromfile.txt yeast-pg/chrom-alignments --alignOptions "--pangenome --reference S288C --outVG " 
```

The results are a HAL and VG file, along with a `cactus-align` log, for each chromosome:
```
ls -hs yeast-pg/chrom-alignments/
total 67M
 844K chrI.hal        1012K chrIII.vg       1.6M chrIX.hal       8.0K chrVI.hal.log     2.2M chrVII.vg      4.0K chrXI.hal.log     2.9M chrXII.vg        3.5M chrXV.hal
 8.0K chrI.hal.log     1.9M chrII.vg        8.0K chrIX.hal.log   3.0M chrVII.hal        948K chrVI.vg       3.3M chrXII.hal        1.3M chrXI.vg         8.0K chrXV.hal.log
 2.6M chrII.hal        748K chrI.vg         1.4M chrIX.vg        4.0K chrVII.hal.log    1.5M chrV.vg        8.0K chrXII.hal.log    2.6M chrXIV.hal       3.0M chrXVI.hal
 8.0K chrII.hal.log    4.9M chrIV.hal       1.9M chrV.hal        1.6M chrVIII.hal       2.1M chrX.hal       2.5M chrXIII.hal       8.0K chrXIV.hal.log   8.0K chrXVI.hal.log
 1.2M chrIII.hal       8.0K chrIV.hal.log   8.0K chrV.hal.log    4.0K chrVIII.hal.log   4.0K chrX.hal.log   4.0K chrXIII.hal.log   2.2M chrXIV.vg        2.3M chrXVI.vg
 8.0K chrIII.hal.log   4.0M chrIV.vg        1.1M chrVI.hal       1.2M chrVIII.vg        1.9M chrXI.hal      1.7M chrXIII.vg        1.7M chrX.vg          2.7M chrXV.vg
```

### Yeast: Joining the Chromosome Alignments

As in the primates example, `cactus-graphmap-join` is used to make the final indexes.  Its use is identical, except multiple graphs are passed as input.  We also pass in the HAL files so it can merge them too.

```
cactus-graphmap-join ./jobstore --vg yeast-pg/chrom-alignments/*.vg --hal yeast-pg/chrom-alignments/*.hal --outDir ./yeast-pg --outName yeast-pg --reference S288C --vcf --giraffe clip

```

The GFA, VCF and all `vg giraffe` indexes will now be in `yeast-pg`:
```
ls -hs yeast-pg/yeast-pg*
 21M yeast-pg/yeast-pg.dist      103M yeast-pg/yeast-pg.min             4.3M yeast-pg/yeast-pg.vcf.gz
 38M yeast-pg/yeast-pg.full.hal  4.6M yeast-pg/yeast-pg.raw.vcf.gz      8.0K yeast-pg/yeast-pg.vcf.gz.tbi
 16M yeast-pg/yeast-pg.gbz       8.0K yeast-pg/yeast-pg.raw.vcf.gz.tbi
 14M yeast-pg/yeast-pg.gfa.gz     16K yeast-pg/yeast-pg.stats.tgz

zcat yeast-pg/yeast-pg.gfa.gz | grep '^W' | awk '{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 }' | grep S288C
W	S288C	0	chrIII	0	341580
W	S288C	0	chrII	0	813597
W	S288C	0	chrI	0	219929
W	S288C	0	chrIV	0	1566853
W	S288C	0	chrIX	0	440036
W	S288C	0	chrVIII	0	581049
W	S288C	0	chrVII	0	1091538
W	S288C	0	chrVI	0	271539
W	S288C	0	chrV	0	583092
W	S288C	0	chrXIII	0	930506
W	S288C	0	chrXII	0	1075542
W	S288C	0	chrXI	0	666862
W	S288C	0	chrXIV	0	777615
W	S288C	0	chrX	0	751611
W	S288C	0	chrXVI	0	954457
W	S288C	0	chrXV	0	1091343

halStats yeast-pg/yeast-pg.full.hal --sequenceStats S288C
SequenceName, Length, NumTopSegments, NumBottomSegments
chrI, 219929, 5766, 0
chrII, 813597, 15896, 0
chrIII, 341580, 7778, 0
chrIV, 1566853, 32488, 0
chrIX, 440036, 10368, 0
chrV, 583092, 12226, 0
chrVI, 271539, 6587, 0
chrVII, 1091538, 19912, 0
chrVIII, 581049, 9996, 0
chrX, 751611, 12635, 0
chrXI, 666862, 11861, 0
chrXII, 1075542, 21991, 0
chrXIII, 930506, 16183, 0
chrXIV, 777615, 17052, 0
chrXV, 1091343, 23011, 0
chrXVI, 954457, 19805, 0
```

### Yeast: Making a UCSC Genome Browser Assembly Hub

The HAL file can be used to produce an assembly hub has follows.  Note that `PYTHONPATH` must set as described in Cactus's installation instructions. 

```
hal2assemblyHub.py ./jobstore ./yeast-pg/yeast-pg.full.hal yeast-pg/hub --shortLabel yeast --longLabel "yeast pangenome"
```

Move `yeast-pg/hub` to somewhere web-accessible, and pass the full URL of `yeast-pg/hub/hub.txt` to the Genome Browser in the "My Data -> Track Hubs" menu.   Select `S288C` as the reference and display the hub.  Right-click on the display and select "Configure yeast track set" to toggle on all the assemblies (and toggle off Anc0 and _MINIGRAPH_).

## HPRC Graph
The [Human Pangenome Reference Consortium](https://humanpangenome.org/data-and-resources/) is producing an ever-growing number of high quality phased assemblies.  This section will demonstrate how to use the Minigraph-Cactus Pangenome Pipeline to construct a Pangenome from them.  Note the instructions here are slightly different than were used to create the v1.0 Minigraph-Cactus pangenome that's been released by the HPRC, as they are based on a more recent and improved version of the pipeline. 

The steps below are run on AWS/S3, and assume everything is written to s3://MYBUCKET. All jobs are run on r5.8xlarge (32 cores / 256G RAM) nodes. In theory, the entire pipeline could therefore be run on a single machine (ideally with 64 cores).  It would take several days though. They can be run on other batch systems, at least in theory.  Most of the compute-heavy tasks spawn relatively few jobs, and may be amenable to SLURM environments.

The following environment variables must be defined: `MYBUCKET` and `MYJOBSTORE`. All output will be placed in `MYBUCKET`, and `MYJOBSTORE` will be used by TOIL for temporary storage.  For example
```
export MYBUCKET=s3://vg-k8s/vgamb/wg/cactus/GRCh38-f1g-90/may4
export MYJOBSTORE=aws:us-west-2:cactus-hprc-jobstore
export VERSION=may4
export MINIGRAPH=https://zenodo.org/record/6499594/files/GRCh38-90c.r518.gfa.gz
```

WDL / cactus-prepare support is in progress!

### HPRC Graph: Setup and Name Munging

**Important** The Cactus-Minigraph Pipeline does not support alt contigs in the reference.  If you really want them in your graph, then you will need to pull them out into separate samples (ie one alt contig per region per sample).  Otherwise they will end up as separate reference contigs and not align together.  As such we advice using the GRCh38 fasta file referenced in the `hprc-${VERSION}-mc.seqfile` generated below for any graph using GRCh38 as a reference.  

The fasta sequences for the Year-1 HPRC assemblies are [available here](https://github.com/human-pangenomics/HPP_Year1_Assemblies).  We begin by using them to create an input seqfile for Cactus:

```
wget -q https://raw.githubusercontent.com/human-pangenomics/HPP_Year1_Assemblies/main/assembly_index/Year1_assemblies_v2_genbank.index
grep GRCh38 Year1_assemblies_v2_genbank.index | sed -e 's/_no_alt_analysis_set\t/\t/g' | awk '{print $1 "\t" $2}' > hprc-${VERSION}-mc.seqfile
printf "CHM13v2\thttps://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/analysis_set/chm13v2.0_maskedY_rCRS.fa.gz\n" >> hprc-${VERSION}-mc.seqfile
tail -n +2 Year1_assemblies_v2_genbank.index | awk '{print $1 ".1\t" $2}' | grep -v CHM13 | grep -v GRCh38 >> hprc-${VERSION}-mc.seqfile
tail -n +2 Year1_assemblies_v2_genbank.index | awk '{print $1 ".2\t" $3}' | grep -v CHM13 | grep -v GRCh38 >> hprc-${VERSION}-mc.seqfile
sort -k1 hprc-${VERSION}-mc.seqfile > hprc-${VERSION}-mc.seqfile.sort ; mv hprc-${VERSION}-mc.seqfile.sort hprc-${VERSION}-mc.seqfile
sed hprc-${VERSION}-mc.seqfile -i -e 's%s3://human-pangenomics/working/%https://s3-us-west-2.amazonaws.com/human-pangenomics/working/%g'
```

We have been holding out three samples for evaluation.  This is, of course, optional, but can be done here with:
```
grep -v 'HG002\|HG005\|NA19240' hprc-${VERSION}-mc.seqfile > t && mv t hprc-${VERSION}-mc.seqfile
```

Also, a misjoin in `HG02080#1#JAHEOW010000073.1` was manually corrected by using `samtools faidx` to break it into `HG02080#1#JAHEOW010000073.1_sub_0_7238466` and `HG02080#1#JAHEOW010000073.1_sub_7238466_12869124`.  The `sub_X_Y` (0-based, open-ended like BED) coordinates are understood by the pipeline, and the offsets will be preserved in the GFA W-lines at the end.  If we don't apply this change, then path names with ":"'s will end up in the HAL which will prevent it from working with assembly hubs.  
```
wget -q $(grep HG02080\.1 hprc-${VERSION}-mc.seqfile | tail -1 | awk '{print $2}') -O HG02080.1.fa.gz
gzip -d HG02080.1.fa.gz
samtools faidx HG02080.1.fa
keep_contigs=$(awk '{print $1}' HG02080.1.fa.fai | grep -v JAHEOW010000073\.1)
samtools faidx HG02080.1.fa ${keep_contigs} > HG02080.1.fix.fa
samtools faidx HG02080.1.fa "HG02080#1#JAHEOW010000073.1:1-7238466" | sed -e 's/\([^:]*\):\([0-9]*\)-\([0-9]*\)/echo "\1_sub_$((\2-1))_\3"/e' >> HG02080.1.fix.fa
samtools faidx HG02080.1.fa "HG02080#1#JAHEOW010000073.1:7238467-12869124" | sed -e 's/\([^:]*\):\([0-9]*\)-\([0-9]*\)/echo "\1_sub_$((\2-1))_\3"/e' >> HG02080.1.fix.fa
bgzip HG02080.1.fix.fa --threads 8
aws s3 cp HG02080.1.fix.fa.gz ${MYBUCKET}/fasta/
grep -v HG02080\.1 hprc-${VERSION}-mc.seqfile > t && mv t hprc-${VERSION}-mc.seqfile
printf "HG02080.1\t${MYBUCKET}/fasta/HG02080.1.fix.fa.gz\n" >> hprc-${VERSION}-mc.seqfile
```

```
head -4 hprc-${VERSION}-mc.seqfile
CHM13v2 https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/analysis_set/chm13v2.0_maskedY_rCRS.fa.gz
GRCh38  https://s3-us-west-2.amazonaws.com/human-pangenomics/working/HPRC_PLUS/GRCh38/assemblies/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz
HG00438.1       https://s3-us-west-2.amazonaws.com/human-pangenomics/working/HPRC/HG00438/assemblies/year1_f1_assembly_v2_genbank/HG00438.paternal.f1_assembly_v2_genbank.fa.gz
HG00438.2       https://s3-us-west-2.amazonaws.com/human-pangenomics/working/HPRC/HG00438/assemblies/year1_f1_assembly_v2_genbank/HG00438.maternal.f1_assembly_v2_genbank.fa.gz
```

The names in these fasta files are for the form `chr1, chr2, etc` in CHM13 and GRCh38, and 

```
SAMPLE#HAPLTOYPE#CONTIG
```

in the other samples.  The "#" symbols cannot be displayed in the UCSC Genome Browser, so it is recommended to stick to the conventions described above: where the fasta contig names are just the `CONTIG`, and the genome name is `SAMPLE.HAPLOTYPE`. Sequence names of the form `SAMPLE#HAPLOTYPE#CONTIG` will be replaced by default with `id=GENOME|CONTIG` by default by `cactus-preprocess --pangenome`.  

We first setup a place for the renamed fasta files using `cactus-prepare` to generate a new seqfile, hprc-pg/hprc-${VERSION}-mc.seqfile
```
cactus-prepare ./hprc-${VERSION}-mc.seqfile --outDir hprc-pg --seqFileOnly
# when running on AWS, data needs to be in S3
sed hprc-pg/hprc-${VERSION}-mc.seqfile -e "s%hprc-pg%${MYBUCKET}/fasta%g" | grep -v ';' > hprc-${VERSION}-mc.pp.seqfile

# save them
aws s3 cp hprc-${VERSION}-mc.seqfile ${MYBUCKET}/
aws s3 cp hprc-${VERSION}-mc.pp.seqfile ${MYBUCKET}/

# finally, we run cactus-preprocess --pangenome
cactus-preprocess ${MYJOBSTORE} hprc-${VERSION}-mc.seqfile hprc-${VERSION}-mc.pp.seqfile --pangenome  --batchSystem mesos --provisioner aws --defaultPreemptable --nodeType r5.8xlarge --nodeStorage 500 --maxNodes 2 --logFile hprc-${VERSION}-mc-grch38.pp.log

```

### HPRC Graph: Mapping to the Graph

Note: since there is already a minigraph available for this data, we just use it instead of constructing it ourselves. See the previous examples for how to construct a minigraph with `cactus-minigraph`.

Now that the sequences are ready, we run `cactus-graphmap` as before.  There is a new option:

`--delFilter N` : Filter out mappings that would induce a deletion bubble of `>N` bases w.r.t. a path in the reference.  If this option is used, the unfiltered paf will also be output (with a `.unfiltered` suffix) as well as a log detailing what was filtered and why (`.filter.log` suffix).  This option is very important as minigraph will produce a small number of split-mappings that can cause chromosome-scale bubbles.  By default, it is set to 1000000.

```
cactus-graphmap ${MYJOBSTORE} hprc-${VERSION}-mc.pp.seqfile ${MINIGRAPH} ${MYBUCKET}/hprc-${VERSION}-mc-grch38.paf --outputGAFDir ${MYBUCKET}/gaf-hprc-${VERSION}-mc-grch38 --outputFasta ${MYBUCKET}/fasta/minigraph.grch38.sv.gfa.fa.gz --reference GRCh38 --mapCores 16 --delFilter 10000000  --batchSystem mesos --provisioner aws --defaultPreemptable --nodeType r5.8xlarge:1.5 --nodeStorage 650 --maxNodes 25 --betaInertia 0 --targetTime 1  --logFile hprc-${VERSION}-mc-grch38.paf.log
```

Note:  The `--betaInertia 0 --targetTime 1` options force Toil to create AWS instances as soon as they are needed.

This command uses the spot market by specifying `:1.35` after the node type to bid $1.35/hr (on-demand pricing at time of writing is about $2.00).

### HPRC Graph: Splitting by Chromosome

There are too many reference contigs to make a graph for each because of all the unplaced contigs in GRCh38.  Ideally, we would drop them but it simplifies some downstream pipelines that use tools that expect them to be in BAM headers etc. to just include them in the graph.  To do this, we use the `--otherContig` option to lump them all into a single job, and `--refContigs` to spell out all the contigs we want to treat separately.  Note that the final output will be the same whether or not `--otherContig` is used. This option serves only to reduce the number of output files (and therefore alignment jobs). 

```
cactus-graphmap-split ${MYJOBSTORE}  hprc-${VERSION}-mc.pp.seqfile ${MINIGRAPH} ${MYBUCKET}/hprc-${VERSION}-mc-grch38.paf --outDir ${MYBUCKET}/chroms-hprc-${VERSION}-mc-grch38 --otherContig chrOther --refContigs $(for i in `seq 22`; do echo chr$i; done ; echo "chrX chrY chrM") --reference GRCh38  --batchSystem mesos --provisioner aws --defaultPreemptable --nodeType r5.8xlarge --nodeStorage 1000 --maxNodes 5 --betaInertia 0 --targetTime 1 --logFile hprc-${VERSION}-mc-grch38.split.log
```

### HPRC Graph: Batch Alignment

The rest of the pipeline is proceeds as in the yeast example. We need to manually download the chromfile though.  We also use a new option

`--maxLen N` : Do not attempt to align more than `N` bases with the Cactus base aligner (activated with `--base`).  This will save aligning too far into anchorless regions, which cannot be properly resolved with base alignment alone.  It is 1000000 by default. 

This command will create a vg and hal file for each chromosome in ${MYBUCKET}/align-batch-grch38/
```

aws s3 cp ${MYBUCKET}/chroms-hprc-${VERSION}-mc-grch38/chromfile.txt .
cactus-align-batch ${MYJOBSTORE} ./chromfile.txt ${MYBUCKET}/align-hprc-${VERSION}-mc-grch38 --alignCores 16  --alignOptions "--pangenome --maxLen 10000 --reference GRCh38   --outVG" --batchSystem mesos --provisioner aws --defaultPreemptable --nodeType r5.8xlarge:1.5 --nodeStorage 1000 --maxNodes 20 --betaInertia 0 --targetTime 1 --logFile hprc-${VERSION}-mc-grch38.align.log
```

### HPRC Graph: Creating the Whole-Genome Graph

Important:
* we use `--filter 9` and `--giraffe` to make giraffe indexes on the subgraph covered by at least 9/90 haplotypes (filter does not apply to the reference).
* we use `--reference GRCh38 CHM13v2` to specify an additional reference. In this case `CHM13v2` will still be clipped (but not filtered), and it will be treated as a reference path in vg (and therefore easier to query).  You can add `--vcfReference GRCh38 CHM13v2` to make a VCF based on CHM13 too.
* we use `--indexCores 63` to allow indexing and hal merging to be done in parallel, which can save quite a bit of time. 

```
cactus-graphmap-join ${MYJOBSTORE} --vg $(for j in $(for i in `seq 22`; do echo chr$i; done ; echo "chrX chrY chrM chrOther"); do echo ${MYBUCKET}/align-hprc-${VERSION}-mc-grch38/${j}.vg; done) --hal $(for j in $(for i in `seq 22`; do echo chr$i; done ; echo "chrX chrY chrM chrOther"); do echo ${MYBUCKET}/align-hprc-${VERSION}-mc-grch38/${j}.hal; done) --outDir ${MYBUCKET}/ --outName hprc-${VERSION}-mc-grch38 --reference GRCh38 --filter 9 --giraffe --vcf --gbz --gfa --vg-chroms --batchSystem mesos --provisioner aws --defaultPreemptable --nodeType r5.16xlarge --nodeStorage 1000 --maxNodes 1 --indexCores 63  --logFile hprc-${VERSION}-mc-grch38.join.log 
```

**All sequences clipped out by `cactus-graphmap-join` will be saved in a BED file in the ".stats.gz" file in its output directory.**

### HPRC Graph: Changing the Reference

The selection of the reference genome is very important, as it will be used as the backbone for the graph.  It is the only genome that is guaranteed to not have any cycles nor to ever be clipped, and therefore provides a coordinate system in the graph.  Any input genome can be used as a reference, provided it's consistently passed as the `--reference` option to all the commands.  It also must not have a "." in its genome name.  In practice, there are usually two possible references for the HPRC graphs: GRCh38 and CHM13. 

It is advisable to also pass `--reference CHM13v2 GRCh38 --vcfReference CHM13v2 GRCh38` to `cactus-graphmap-join` to tell it to make a second VCF based on the GRCh38 reference.

Note: some contig names like `chrY` (if it is not included) and options like `--otherContig` will not be necessary for `CHM13`

### HPRC Graph: Other Approaches for Masking or Clipping out Complex Regions

The Pangenome Pipeline supports options to for special handling of masked regions at pretty much every step.  These were added to address various issues during initial development and testing.  The approach of just aligning everything and filtering based on the minigraph described above is much simpler and seems at least as effective.  

**cactus-preprocess**

Most satellite sequence can be detected with `dna-brnn`, which can be run with via the `--maskAlpha --minLength 100000 --brnnCores 8` options in `cactus-preprocess`.  The entire pipeline supports sub-sequence fragments via naming conventions, so the masked sequence can be clipped out instead of masked by using `--clipAlpha` instead of `--maskAlpha`

The minigraph mappings themselves can also be used to derive regions to mask, by finding gaps in the alignments.  This can be done by passing a PAF file (output from `cactus-graphmap`) back into `cactus-preprocess` via the `--maskFile` option.  This option can also accept BED files to mask any user-specified regions.  When using this option, the `--maskAction` option can be used to specify whether masked sequence is clipped out or not.

**cactus-graphmap**

Softmasked input sequence can be ignored by using the `--maskFilter 100000` option.  This will force such sequence to remain unaligned.

**cactus-graphmap-split**

Softmasked input sequence can (and should) be ignored when computing coverage in order to assign contigs to reference chromosomes.  This is done with `--maskFilter 100000`

**cactus-align-batch**

Softmasked input can be ignored (and forced to stay unaligned) with the `--barMaskFilter 100000` option to `cactus-align`, or by including it in the `cactus-align-batch --alignOptions "--barMaskFilter 100000"`

### HPRC Version 1.0 Graphs

These graphs were created with the [cactus-pangenome.sh script](https://github.com/glennhickey/pg-stuff/blob/c87b9236a20272b127ea2fadffc5428c5bf15c0e/cactus-pangenome.sh) using Cactus commit [6cd9a42cdf40ad61843664ed82c9d5bc26445570](https://github.com/ComparativeGenomicsToolkit/cactus/commit/6cd9a42cdf40ad61843664ed82c9d5bc26445570).  The seqfile input was constructed as above, except chrY was only added to CHM13 for the CHM13-based graph (and chrEBV was never added).  Instead, a decoy graph consisting of chrEBV and all the hs38d1 contigs was added to both graphs in the `cactus-graphmap-join` step.

The other main differences between this pipeline and ${VERSION} are
* Input fasta files were softmasked with dna-brnn regions `>100kb`
* After mapping to the graph, minimizer gaps `>100kb` were masked using a second call to `cactus-preprocess`
* The two sets of masked regions were merged together and clipped out the input sequences.
* The clipped sequences were remapped to the graph once again and the pipeline continued from there
* The `--base` option was never used to perform sequence-to-graph base alignment (it didn't exist)
* The `--delFilter` option didn't exist either, so several large spurious bubbles made it into the graphs
* Much more stringent options were used to assign contigs to chromosomes with `cactus-graphmap-join`.  This was possible to some extent because the contigs were clipped, but also caused more sequence to be classified as ambiguous.
* `cactus-graphmap-join` clipped out sequence that was unaligned to anything else (including minigraph), rather than unaligned to minigraph.  (this is less stringent).
* A few small bugs were fixed in Cactus between the two versions, notably one that caused erroneous tiny duplications and inversions. 

GRCh38 graph command line
```
./cactus-pangenome.sh -j aws:us-west-2:glennhickey-jobstore7 -s ./hprc-year1-f1g.fix.HG02080.1.brnn.leaveout.seqfile -m ftp://ftp.dfci.harvard.edu/pub/hli/minigraph/HPRC-f1g/GRCh38-f1g-90.gfa.gz  -o s3://vg-k8s/vgamb/wg/cactus/GRCh38-f1g-90/aug11 -n GRCh38-f1g-90-mc-aug11 -r GRCh38 -d s3://vg-k8s/vgamb/wg/fasta/hs38d1.decoys.only.vg  -g  -F  -C -M 100000 -K 10000  2>> stderr.aug11.2.log > /dev/null
```

CHM13 graph command line
```
 ./cactus-pangenome.sh -j aws:us-west-2:glennhickey-jobstore-hprc4 -s ./hprc-year1-f1g.chmy.fix.HG02080.1.brnn.leaveout.seqfile -m ftp://ftp.dfci.harvard.edu/pub/hli/minigraph/HPRC-f1g/CHM13-f1g-90.gfa.gz  -o s3://vg-k8s/vgamb/wg/cactus/CHM13-f1g-90/aug11 -n CHM13-f1g-90-mc-aug11  -r CHM13 -v GRCh38 -d s3://vg-k8s/vgamb/wg/fasta/hs38d1.decoys.only.vg  -g  -F -C -M 100000 -K 10000 -y  2>> stderr.aug11.chm13.3.log > /dev/null
```

## Frequently Asked Questions

Q; `cactus-graphmap-join` keeps crashing with a segfault and/or running out of memory.

A: This is usually because `vg index -j` crashes while computing the distance index for `--giraffe`.  Make sure you did not disable clipping or filtering.  Try to use `--filter N` where `N` is about 10% of your input samples.  If that fails, try on a system with more memory.

Q: Why are the node id's different in my GFA and GBZ

A: GBZ construction chops nodes to a maximum length of 1024, which changes their ids.  The GFA and VCF use the original, unchopped IDs.  You can get a mapping between the chopped and unchopped IDs from the GBZ using `vg gbwt --translation`.  This is all really annoying and I wonder if it's not just better to chop everything?

Q: My tools can't read GFA 1.1 and W-lines.

A: As mentioned [above](#output), you can use `vg convert -gfW` to convert from GFA 1.1 to GFA 1.0.

Q: I hate the idea of clipping and filtering sequence from my graph.  Why do I have to do it?!

A: So current toolchains can work with your graphs.  But clipping and filtering is optional in `cactus-graphmap-join`, just specify `full` for the various outputs (keeping in mind that much satellite sequence won't be aligned to anything).


