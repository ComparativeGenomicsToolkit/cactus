The Minigraph-Cactus Pangenome Pipeline
===

## Introduction

[Cactus](../README.md) uses a phylogenetic tree as a guide in order to progressively create multiple alignments.  This heuristic allows Cactus to scale linearly with the number of input genomes, by decomposing the work into one alignment per internal (ancestral) node of the tree.  If the guide tree is fully resolved (binary), only two genomes (plus up to three outgroups) are aligned in each subproblem.

Progressively aligning up a guide tree makes sense when the evolution of the input genomes can be explained by a tree.  It is robust to small errors in the tree, as well as small numbers of non-treelike events (ex incomplete lineage sorting or horizontal genen transfer), making it [suitable for alignments of different vertebrate species](https://doi.org/10.1038/s41586-020-2871-y).

But the tree-like assumption breaks down when considering an alignment of individuals from the *same species*.  Such within-population genome alignments are increasingly in demand as high-quality assemblies become more available (ex: [HPP](https://humanpangenome.org/)), given their potential to better identify and represent structural variation than more traditional reference-based re-sequencing approaches.

The Minigraph-Cactus Pangenome Pipeline adapts [Cactus](../README.md) to no longer rely on a guide tree, by taking advantage of the relative similarity of the input sequence to use minimizer sketches to determine initial anchors, then partial order alignments to refine them.  It also provides the options to generate output in standard pangenome graph formats such as [vg](https://github.com/vgteam/vg) and [GFA](https://github.com/GFA-spec/GFA-spec), in addition to the usual HAL. 

The other key difference between the Pangenome and Progressive modes is that in Pangenome mode, Cactus will leave one "reference" genome uncollapsed (ie it will never be aligned to itself).  This lets it be used as a unique coordinate system for downstream analyses.  For example, when creating a human pangenome, GRCh38 (or CHM13) can be included and flagged as the reference.  These reference coordinates simplify projection to VCF, allow [vg giraffe](https://doi.org/10.1126/science.abg8871) to produce mappings in BAM format, are compatible with rGFA, etc.

*This is a work in progress and is not yet published.* 

## Overview

The pangenome interface is very similar to the [default utilization of Cactus](../README.md), and is dependent on a [seqfile mapping genome names to fasta locations](seqFile-the-input-file).  The main difference is that a tree need not be provided at the top of the seqfile.  If a tree is present, it must be a star tree (all leaves connected to one root).

If you compiled Cactus yourself as opposed to using a binary or docker release, it is important to run `build-tools/downloadPangenomeTools` to install the pangeonme-specific binaries.

Before running the Minigraph-Cactus Pangenome Pipeline on your data, it is highly recommended that you first try the Evolver Primates and Yeast examples below -- the data is included in Cactus.  These examples will explain some important options and naming conventions.  

### Evolver Primates: Preprocessing and Naming

For genome names in the input seqfile, it is best to use the convention of `SAMPLE.HAPLOTYPE` if possible, where `HAPLOTYPE` is an integer.  Examples of naming within the Seqfile:


```
# Diploid sample:
HG002.1  ./HG002.paternal.fa
HG002.2  ./HG002.maternal.fa

# Haploid sample:
CHM13.0  ./chm13.fa

# Reference Genome (do not specify haplotype):
GRCh38   ./grch38.fa
```

**IMPORTANT:** If genomes do not follow this convention, the paths will not be properly stored in the `vg giraffe` indexes create by `cactus-graphmap-join --giraffe`.

So to begin, we add ".0" to the genome names in the Evolver Primates example and strip out the tree.  We leave "simChimp" alone as it will be our reference:

```
tail -n +2 examples/evolverPrimates.txt | sed -e 's/^simHuman/simHuman.0/g' -e 's/^simGorilla/simGorilla.0/g' -e 's/^simOrang/simOrang.0/g' > evolverPrimates.pg.txt
```

Input contig names must be *globally* unique to avoid collisions in PAF files.  This can be enforced with `cactus-preprocess` by adding a prefix to each one. Ex:

```
# use cactus-prepare to make an output seqfile 
cactus-prepare ./evolverPrimates.pg.txt --outDir primates-pg --seqFileOnly
# fix the contig names
cactus-preprocess ./jobstore ./evolverPrimates.pg.txt primates-pg/evolverPrimates.pg.txt --skipMasking 
```

Note: The pangenome pipeline supports gzipped fasta files if they end with `.gz`, likewise for `.gfa`.  (`.paf` files must be left uncompressed)

### Evolver Primates: Constructing the Minigraph GFA

The next step is to create the initial [minigraph](https://github.com/lh3/minigraph) graph. `minigraph` iteratively adds sequences to the graph, beginning with a reference genome that will serve throughout the pipeline as its uncollapsed backbone.  This reference genome needs to selected now and used consistently (via the `--reference` option in all following commands).  Other genomes are added to the `minigraph` in the order they appear in the seqfile.

```
cactus-minigraph ./jobstore primates-pg/evolverPrimates.pg.txt primates-pg/primates.gfa.gz --realTimeLogging --reference simChimp
```

### Evolver Primates: Mapping the Genomes Back to the Minigraph

`minigraph` does not perform base alignment, so the graph constructed above will only contain SVs (>50bp).  The remainder of the pipeline will add base-level alignments to the graph using Cactus.  The first step is to use `minigraph` to map each input sequence back to the graph.  The `minigraph` GFA itself will become a "genome" in the Cactus graph, so a path for its fasta sequence also needs to be specified with `--outputFasta`.  

```
cactus-graphmap ./jobstore primates-pg/evolverPrimates.pg.txt primates-pg/primates.gfa.gz primates-pg/primates.paf --realTimeLogging --reference simChimp --outputFasta primates-pg/primates.gfa.fa.gz
```

You can tune the number of CPUs used by each mapping job with the `--mapCores` option.  This command produces a mapping of input sequences to the graph sequences in PAF format with cigars.  

### Evolver Primates: Creating the Cactus Alignment

Now Cactus can be run, much as it would its default progressive mode.  The only difference is that we specify the `--pangenome`, `--outVG` and `--reference` options.

```
cactus-align ./jobstore primates-pg/evolverPrimates.pg.txt primates-pg/primates.paf primates-pg/primates.hal --pangenome --outVG --reference simChimp --realTimeLogging
```

### Evolver Primates: Creating the VG Indexes

Finally, indexes for `vg giraffe` as well as a VCF file can be created.  **Important:** For the paths to be properly represented, the genomes need to be named in `SAMPLE.HAPLOTYPE` format as described above, and `--wlineSep .` must be passed.

The `--gfaffix` option zips up redundant bubbles and is also important to use.

```
cactus-graphmap-join ./jobstore --vg primates-pg/primates.vg --outDir ./primates-pg --outName primates-pg --reference simChimp --vcf --giraffe --gfaffix  --wlineSep "." --realTimeLogging
```

If it worked properly, the reference contigs, prefixed by the genome name and a "." will be in the xg, as well as P-lines in the GFA:
```
vg paths -Ev primates-pg/primates-pg.xg
simChimp.simChimp.chr6	596350

zcat primates-pg/primates-pg.gfa.gz | grep ^P | awk '{print $1 "\t" $2}'
P	simChimp.simChimp.chr6
```

and all other paths will be listed as haplotype threads.  They will be stored in the GBWT, as well as in W-lines in the GFA:
```
vg paths -Lg  primates-pg/primates-pg.gbwt
_thread__gbwt_ref_simChimp.simChimp.chr6_0_0
_thread_simGorilla_simGorilla.chr6_0_0
_thread_simHuman_simHuman.chr6_0_0
_thread_simOrang_simOrang.chr6_0_0

zcat primates-pg/primates-pg.gfa.gz | grep ^W | awk '{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6}' 
W	simGorilla	0	simGorilla.chr6	0	599081
W	simHuman	0	simHuman.chr6	0	597871
W	simOrang	0	simOrang.chr6	0	591073
```

The VCF will be based on the reference path (simChimp) and have a sample for each haplotype :
```
gzip -dc primates-pg/primates-pg.vcf.gz | grep CHROM -A 1
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	simGorilla	simHuman	simOrang
simChimp.simChimp.chr6	15	>1>4	T	G	60	.	AC=1;AF=0.333333;AN=3;AT=>1>2>4,>1>3>4;NS=3;LV=0	GT	0	0	1
```

Reference paths are needed for coordinates, but are inefficient to store in vg.  Haplotype threads are exploited by `vg giraffe` and can be efficiently compressed in the GBWT.

## Splitting by Chromosome

`cactus-align`'s memory usage will become prohibitive after about `10G` of input sequence. So in order to scale to dozens of human-sized genomes, `cactus-align` must be run individually on each chromosome.  By definition, inter-chromosomal events will *not be* represented using this approach.

### Yeast: Getting Started

Below is an example of creating a yeast pangenome chromosome by chromosome, referenced on S288C.  The genome names are already in SAMPLE.HAPLOTYPE format, and the first steps are identical to above.  

```
# make the seqfile
cactus-prepare ./examples/yeastPangenome.txt --outDir yeast-pg --seqFileOnly

# add unique prefixes to fasta contigs
cactus-preprocess ./jobstore ./examples/yeastPangenome.txt ./yeast-pg/yeastPangenome.txt  --skipMasking

# make the minigraph
cactus-minigraph ./jobstore  ./yeast-pg/yeastPangenome.txt ./yeast-pg/yeast.gfa --realTimeLogging --reference S288C

# map back to the minigraph
cactus-graphmap ./jobstore ./yeast-pg/yeastPangenome.txt ./yeast-pg/yeast.gfa ./yeast-pg/yeast.paf --outputFasta ./yeast-pg/yeast.gfa.fa --realTimeLogging --reference S288C
```

### Yeast: Splitting By Chromosome
Now the PAF and GFA `minigraph` output can be used to partition the graph and mappings based on the reference genome's (S288C's) chromosomes:

```
cactus-graphmap-split ./jobstore ./yeast-pg/yeastPangenome.txt ./yeast-pg/yeast.gfa ./yeast-pg/yeast.paf --outDir yeast-pg/chroms --realTimeLogging --reference S288C
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
zcat yeast-pg/chroms/_AMBIGUOUS_/fasta/UWOPS034614.0__AMBIGUOUS_.fa.gz | grep '>'
>id=UWOPS034614.0|chrXIII
>id=UWOPS034614.0|chrX
>id=UWOPS034614.0|chrVII
>id=UWOPS034614.0|chrXI
>id=UWOPS034614.0|chrVIII
```

The reason why these contigs are unassigned to a chromosome can normally be found in `minigraph.split.log`:

```
grep ambiguous yeast-pg/chroms/minigraph.split.log -A 3
Query contig is ambiguous: id=UWOPS034614.0|chrXI  len=792116 cov=0.551694 (vs 0.5) uf=1.41798 (vs 3)
 Reference contig mappings:
  chrVII: 308189
  chrXI: 437006
--
Query contig is ambiguous: id=UWOPS034614.0|chrVIII  len=738767 cov=0.47896 (vs 0.5) uf=1.04615 (vs 3)
 Reference contig mappings:
  chrVII: 353840
  chrVIII: 338232
--
Query contig is ambiguous: id=UWOPS034614.0|chrXIII  len=662343 cov=0.526851 (vs 0.5) uf=1.72564 (vs 3)
 Reference contig mappings:
  chrVII: 91734
  chrXI: 202218
--
Query contig is ambiguous: id=UWOPS034614.0|chrVII  len=632616 cov=0.409076 (vs 0.5) uf=1.54897 (vs 3)
 Reference contig mappings:
  chrVII: 258788
  chrVIII: 167071
--
Query contig is ambiguous: id=UWOPS034614.0|chrX  len=1092164 cov=0.486244 (vs 0.25) uf=1.05286 (vs 3)
 Reference contig mappings:
  chrX: 504394
  chrXIII: 531058

```

It is saying that not enough bases in these contigs aligned to a single reference chromosome, given the 50% threshold and 3X uniqueness factor.  These thresholds are explained, and can be adjusted in, the `<graphmap-split>` section of the [cactus config](../src/cactus/cactus_progressive_config.xml).

The sequence-to-graph PAF files can be refined by adding the `--base` option to `cactus-graphmap`.  This can often improve the accuracy of `cactus-graphmap-split`.  This option will be used in the HPRC example below. 

### Yeast: Batch Aligning the Chromosomes

`cactus-align` can be run individually on each chromosome using the seqFiles created above.  This can be automated using the `cactus-align-batch` script, which takes as input a "chromFile", which is just a list of seqFiles and PAFs.  Such a chromFile was generated by `cactus-graphmap-split` above.

The options we would normally pass directly to `cactus-align` must be quoted and passed via `--alignOptions` here:

```
cactus-align-batch ./jobstore ./yeast-pg/chroms/chromfile.txt yeast-pg/chrom-alignments --alignOptions "--pangenome --reference S288C --outVG --realTimeLogging" --realTimeLogging
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
cactus-graphmap-join ./jobstore --vg yeast-pg/chrom-alignments/*.vg --hal yeast-pg/chrom-alignments/*.hal --outDir ./yeast-pg --outName yeast-pg --reference S288C --vcf --giraffe --gfaffix  --wlineSep "." --realTimeLogging

```

The GFA, VCF and all `vg giraffe` indexes will now be in `yeast-pg`:
```
ls -hs yeast-pg/yeast-pg*
 20M yeast-pg/yeast-pg.dist     37M yeast-pg/yeast-pg.gg   860K yeast-pg/yeast-pg.snarls    8.0K yeast-pg/yeast-pg.vcf.gz.tbi
7.6M yeast-pg/yeast-pg.gbwt     36M yeast-pg/yeast-pg.hal  3.1M yeast-pg/yeast-pg.trans.gz   24M yeast-pg/yeast-pg.xg
 14M yeast-pg/yeast-pg.gfa.gz  111M yeast-pg/yeast-pg.min  4.8M yeast-pg/yeast-pg.vcf.gz

vg paths -Ev yeast-pg/yeast-pg.xg
S288C.chrIII	341580
S288C.chrII	813597
S288C.chrI	219929
S288C.chrIV	1566853
S288C.chrIX	440036
S288C.chrVIII	581049
S288C.chrVII	1091538
S288C.chrVI	271539
S288C.chrV	583092
S288C.chrXIII	930506
S288C.chrXII	1075542
S288C.chrXI	666862
S288C.chrXIV	777615
S288C.chrX	751611
S288C.chrXVI	954457
S288C.chrXV	1091343

zcat yeast-pg/yeast-pg.gfa.gz | grep ^P | awk '{print $1 "\t" $2}'
P	S288C.chrIII
P	S288C.chrII
P	S288C.chrI
P	S288C.chrIV
P	S288C.chrIX
P	S288C.chrVIII
P	S288C.chrVII
P	S288C.chrVI
P	S288C.chrV
P	S288C.chrXIII
P	S288C.chrXII
P	S288C.chrXI
P	S288C.chrXIV
P	S288C.chrX
P	S288C.chrXVI
P	S288C.chrXV

halStats yeast-pg/yeast-pg.hal --sequenceStats S288C
SequenceName, Length, NumTopSegments, NumBottomSegments
S288C.chrI, 219929, 3533, 0
S288C.chrII, 813597, 13395, 0
S288C.chrIII, 341580, 6207, 0
S288C.chrIV, 1566853, 27433, 0
S288C.chrIX, 440036, 8115, 0
S288C.chrV, 583092, 9739, 0
S288C.chrVI, 271539, 5140, 0
S288C.chrVII, 1091538, 17276, 0
S288C.chrVIII, 581049, 8089, 0
S288C.chrX, 751611, 11173, 0
S288C.chrXI, 666862, 10212, 0
S288C.chrXII, 1075542, 18179, 0
S288C.chrXIII, 930506, 13978, 0
S288C.chrXIV, 777615, 14338, 0
S288C.chrXV, 1091343, 19735, 0
S288C.chrXVI, 954457, 17092, 0

vg paths -Lg  yeast-pg/yeast-pg.gbwt | head -4
_thread__gbwt_ref_S288C.chrIII_0_0
_thread_DBVPG6044_chrIII_0_0
_thread_SK1_chrIII_0_0
_thread_UWOPS034614_chrIII_0_0

zcat yeast-pg/yeast-pg.gfa.gz | grep ^W | awk '{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6}'  | head -4
zcat yeast-pg/yeast-pg.gfa.gz | grep ^W | awk '{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6}'  | head -4
W	DBVPG6044	0	chrIII	0	332771
W	SK1	0	chrIII	0	340914
W	UWOPS034614	0	chrIII	0	309137
W	Y12	0	chrIII	0	322503
```

### Yeast: Making a UCSC Genome Browser Assembly Hub

The HAL file can be used to produce an assembly hub has follows.  Note that `PYTHONPATH` must set as described in Cactus's installation instructions. 

```
hal2assemblyHub.py ./jobstore ./yeast-pg/yeast-pg.hal yeast-pg/hub --shortLabel yeast --longLabel "yeast pangenome"
```

Move `yeast-pg/hub` to somewhere web-accessible, and pass the full URL of `yeast-pg/hub/hub.txt` to the Genome Browser in the "My Data -> Track Hubs" menu.   Select `S288C` as the reference and display the hub.  Right-click on the display and select "Configure yeast track set" to toggle on all the assemblies (and toggle off Anc0 and _MINIGRAPH_).

## HPRC Graph

The [Human Pangenome Reference Consortium](https://humanpangenome.org/data-and-resources/) is producing an ever-growing number of high quality phased assemblies.  This section will demonstrate how to use the Cactus-Minigraph Pangenome Pipeline to construct a Pangenome from them.  Note the instructions here are slightly different than were used to create the v1.0 Cactus-Minigraph pangenome that's been released by the HPRC, as they are based on a more recent and improved version of the pipeline. 

The steps below are run on AWS/S3, and assume everything is written to s3://MYBUCKET. They can be run on other batch systems, at least in theory.  Most of the compute-heavy tasks spawn relatively few jobs, and may be amenable to SLURM environments.  

WDL / cactus-prepare support is in progress!

### HPRC Graph: Setup and Name Munging

The fasta sequences for the Year-1 HPRC assemblies are [available here](https://github.com/human-pangenomics/HPP_Year1_Assemblies).  Of note, the contig names are of the format:

```
SAMPLE#HAPLTOYPE#CONTIG
```

"#" symbols cannot be displayed in the UCSC Genome Browser, so it is recommended to stick to the conventions described above: where the fasta contig names are just the `CONTIG`, and the genome name is `SAMPLE.HAPLOTYPE`.  This can be accomplished with `cactus-preprocess` (note: the full path of `cactus_progressive_config.xml` will be dependent on your cactus installation directory, and may need to be adjusted):

```
sed src/cactus/cactus_progressive_config.xml -e "s/cutBefore=\"\"/cutBefore=\"#\"/g" > config_cut_hash.xml
```

Now the seqFile can be made as follows:
```
wget https://raw.githubusercontent.com/human-pangenomics/HPP_Year1_Assemblies/main/assembly_index/Year1_assemblies_v2_genbank.index
grep GRCh38 Year1_assemblies_v2_genbank.index | sed -e 's/_no_alt_analysis_set\t/\t/g' | awk '{print $1 "\t" $2}' > hprc.seqfile
grep CHM13 Year1_assemblies_v2_genbank.index | sed -e 's/CHM13_v1.1/CHM13_1_1.0/g' | awk '{print $1 "\t" $2}' >> hprc.seqfile
tail -n +2 Year1_assemblies_v2_genbank.index | awk '{print $1 ".1\t" $2}' | grep -v CHM13 | grep -v GRCh38 >> hprc.seqfile
tail -n +2 Year1_assemblies_v2_genbank.index | awk '{print $1 ".2\t" $3}' | grep -v CHM13 | grep -v GRCh38 >> hprc.seqfile
sort -k1 hprc.seqfile > hprc.seqfile.sort ; mv hprc.seqfile.sort hprc.seqfile
sed hprc.seqfile -i -e 's%s3://human-pangenomics/working/%https://s3-us-west-2.amazonaws.com/human-pangenomics/working/%g'
```
Because we are using GRCh38 as a reference, we added the `.0` suffix to CHM13 (and removed the "." from its version number!).

We have been holding out three samples for evaluation.  This is, of course, optional, but can be done here with:
```
grep -v 'HG002\|HG005\|NA19240' hprc.seqfile > t && mv t hprc.seqfile
```

Also, a misjoin in `HG02080#1#JAHEOW010000073.1` was manually corrected by using `samtools faidx` to break it into the following fragments
```
HG02080#1#JAHEOW010000073.1_sub_0_7238466
HG02080#1#JAHEOW010000073.1_sub_7238466_12869124
```
Note: If the `_sub_start_end` naming convention is used for fragments, the coordinates will be properly preserved throughout the entire pipeline.

```
head -4 hprc.seqfile
CHM13_1_1.0	https://s3-us-west-2.amazonaws.com/human-pangenomics/working/HPRC_PLUS/CHM13/assemblies/chm13.draft_v1.1.fasta.gz
GRCh38	https://s3-us-west-2.amazonaws.com/human-pangenomics/working/HPRC_PLUS/GRCh38/assemblies/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz
HG002.1	https://s3-us-west-2.amazonaws.com/human-pangenomics/working/HPRC_PLUS/HG002/assemblies/year1_f1_assembly_v2_genbank/HG002.paternal.f1_assembly_v2_genbank.fa.gz
HG002.2	https://s3-us-west-2.amazonaws.com/human-pangenomics/working/HPRC_PLUS/HG002/assemblies/year1_f1_assembly_v2_genbank/HG002.maternal.f1_assembly_v2_genbank.fa.gz
```

If we forget to do this, genomes can be renamed before the `vg giraffe` indexes are made using the `--rename` option in `cactus-graphmap-join`. This can also be done manually on the hal files with `halRenameGenomes` or the vg files with `clip-vg -r`.

```
cactus-prepare ./hprc.seqfile --outDir hprc-pg --seqFileOnly
# when running on AWS, data needs to be in S3
sed hprc-pg/hprc.seqfile -i -e 's/hprc-pg/s3:\/\/MYBUCKET\/fasta/'

cactus-preprocess aws:us-west-2:MYJOBSTORE hprc.seqfile hprc-pg/hprc.seqfile --configFile ./config_cut_hash.xml --realTimeLogging --batchSystem mesos --provisioner aws --defaultPreemptable --nodeType r5.4xlarge --nodeStorage 500 --maxNodes 1 

```

### HPRC Graph: Mapping to the Graph

Note: since there is already a minigraph available for this data, we just use it instead of constructing it ourselves. See the previous examples for how to construct a minigraph with `cactus-minigraph`.

Now that the sequences are ready, we run `cactus-graphmap` as before.  There are some new options:

`--delFilter N` : Filter out mappings that would induce a deletion bubble of `>N` bases w.r.t. a path in the reference.  If this option is used, the unfiltered paf will also be output (with a `.unfiltered` suffix) as well as a log detailing what was filtered and why (`.filter.log` suffix).  This option is very important as minigraph will produce a small number of split-mappings that can cause chromosome-scale bubbles.

`--base` : `minigraph` only produces minimizer hits, and these can sometimes be misleading when computing coverage for the above filter or during chromosome splitting.  This option sends each minigraph mapping through cactus to fill in base alignments, leading to more accurate PAFs.  It costs more to compute them, but saves time in down the road in `cactus-align` as there will be much fewer bases to align in that step.

`--maxLen N` : Do not attempt to align more than `N` bases with the Cactus base aligner (activated with `--base`).  This will save aligning too far into anchorless regions.

`--mapCores N` : Use `N` cores for each mapping job.

```
cactus-graphmap aws:us-west-2:MYJOBSTORE hprc-pg/hprc.seqfile https://s3-us-west-2.amazonaws.com/human-pangenomics/pangenomes/freeze/freeze1/minigraph/hprc-v1.0-minigraph-grch38.gfa.gz s3://MYBUCKET/hprc.grch38.paf --outputGAFDir s3://MYBUCKET/gaf-hprc.grch38 --outputFasta s3://MYBUCKET/fasta/minigraph.grch38.gfa.fa.gz --reference GRCh38  --delFilter 10000000 --base --maxLen 100000 --mapCores 16 --realTimeLogging --batchSystem mesos --provisioner aws --defaultPreemptable --nodeType r5.8xlarge --nodeStorage 500 --maxNodes 25 --betaInertia 0 --targetTime 1
```

Note:  The `--betaInertia 0 --targetTime 1` options force Toil to create AWS instances as soon as they are needed.  

### HPRC Graph: Splitting by Chromosome

There are too many reference contigs to make a graph for each because of all the unplaced contigs in GRCh38.  Ideally, we would drop them but it simplifies some downstream pipelines that use tools that expect them to be in BAM headers etc. to just include them in the graph.  To do this, we use the `--otherContig` option to lump them all into a single job, and `--refContigs` to spell out all the contigs we want to treat separately.

```
cactus-graphmap-split aws:us-west-2:MYJOBSTORE  hprc-pg/hprc.seqfile https://s3-us-west-2.amazonaws.com/human-pangenomics/pangenomes/freeze/freeze1/minigraph/hprc-v1.0-minigraph-grch38.gfa.gz s3://MYBUCKET/hprc.grch38.paf --outDir s3://MYBUCKET/chroms-hprc.grch38 --otherContig chrOther --refContigs $(for i in `seq 22`; do echo chr$i; done ; echo "chrX chrY chrM") --reference GRCh38 --realTimeLogging --batchSystem mesos --provisioner aws --defaultPreemptable --nodeType r5.8xlarge --nodeStorage 1000 --maxNodes 1
```

### HPRC Graph: Batch Alignment

The rest of the pipeline is proceeds as in the yeast example. We need to manually download the chromfile though.

This command will create a vg and hal file for each chromosome in s3://MYBUCKET/align-batch-grch38/
```
aws s3 cp s3://MYBUCKET/chroms-grch38/chromfile.txt .
cactus-align-batch aws:us-west-2:MYJOBSTORE ./chromfile.txt s3://MYBUCKET/align-batch-grch38 --alignCores 16 --realTimeLogging --alignOptions "--pangenome --reference GRCh38 --realTimeLogging  --outVG --maxLen 100000" --batchSystem mesos --provisioner aws --defaultPreemptable --nodeType r5.8xlarge --nodeStorage 1000 --maxNodes 10 --betaInertia 0 --targetTime 1
```

### HPRC Graph: Creating the Whole-Genome Graph

The individual chromosome graphs can now be merged as follows:

```
cactus-graphmap-join aws:us-west-2:MYJOBSTORE --vg $(for j in $(for i in `seq 22`; do echo chr$i; done ; echo "chrX chrY chrM chrOther"); do echo s3://MYBUCKET/align-batch-grch38/${j}.vg; done) --hal $(for j in $(for i in `seq 22`; do echo chr$i; done ; echo "chrX chrY chrM chrOther"); do echo s3://MYBUCKET/align-batch-grch38/${j}.hal; done) --outDir s3://MYBUCKET/join-grch38 --outName grch38-hprc --reference GRCh38 --gfaffix  --wlineSep "."  --batchSystem mesos --provisioner aws --defaultPreemptable --nodeType r5.8xlarge --nodeStorage 1000 --maxNodes 1 --indexCores 31 --realTimeLogging
```

Note: `--indexCores 31` is used (instead of 32) to let the HAL merging job run in parallel on the same machine. 

### HPRC Graph: Filtering Complex Regions and Indexing for Giraffe

The graph created above will contain all input contigs except those which could not be unambiguously assigned to any chromosome. But there will be a lot of unaligned (and poorly aligned) sequence in and around centromeres, some acrocentric short arms, etc. In most cases, this will lead to long stretches of sequence "dangling" from the graph or forming large bubbles.  These structures are artifacts and should be used cautiously for making inferences about structural variation.  Furthermore, they can render efficient indexing and read mapping impossible with most current tools.

The complex regions can be filtered out with a simple heuristic: remove paths of >Nbp that do not align to the underlying minigraph.  The intuition behind this is to restrict the graph to the SVs present from minigraph (which are generally quite clean and of high confidence), removing anything that arose by failure to map back to the minigraph.  This is done with the `--clipLength N --clipNonMinigraph` options to `graphmap-join`.

We also use the `--giraffe --vcf` options to create the Giraffe indexes and VCF.

We use as input the vg files created by the previous call of `graphmap-join` as well as the `--preserveIDs` option to ensure that our new graph is ID-compatible with the full graph.

```
cactus-graphmap-join aws:us-west-2:MYJOBSTORE --vg $(for j in $(for i in `seq 22`; do echo chr$i; done ; echo "chrX chrY chrM chrOther"); do echo s3://MYBUCKET/clip-grch38-hprc/${j}.vg; done) --outDir s3://MYBUCKET/join-grch38-clip10k --outName grch38-hprc-clip10k --reference GRCh38  --wlineSep "." --clipLength 10000 --clipNonMinigraph --vcf --giraffe --preserveIDs --batchSystem mesos --provisioner aws --defaultPreemptable --nodeType r5.8xlarge --nodeStorage 1000 --maxNodes 1 --indexCores 32 --realTimeLogging
```

**All sequences clipped out by `cactus-graphmap-join` will be saved in BED files in its output directory.**

### HPRC Graph: Filtering by Allele Frequency

It's a work in progress, but the Giraffe-DeepVariant pipeline performs best when further filtering the graph with an allele frequency filter.  Doing so removes rare variants, by definition, but also many assembly and alignment errors.  It can be done using the `--vgClipOpts` option:

```
cactus-graphmap-join aws:us-west-2:MYJOBSTORE --vg $(for j in $(for i in `seq 22`; do echo chr$i; done ; echo "chrX chrY chrM chrOther"); do echo s3://MYBUCKET/clip-grch38-hprc/${j}.vg; done) --outDir s3://MYBUCKET/join-grch38-maf10 --outName grch38-hprc-maf10 --reference GRCh38  --wlineSep "." --vgClipOpts "-d 9 -m 1000" --preserveIDs --batchSystem mesos --provisioner aws --defaultPreemptable --nodeType r5.8xlarge --nodeStorage 1000 --maxNodes 1 --indexCores 32 --realTimeLogging
```

Here `--vgClipOpts "-d 9 -m 1000"` will remove all nodes with fewer than 9 paths covering them, filtering out resulting path fragments of fewer than 1kb bases.

Note: Don't use `--vcf` with this command to make an allele-frequency-filtered VCF.  It is much better to filter the VCF constructed in the previous section directly (ex. with `bcftools`).

### HPRC Graph: Changing the Reference

The selection of the reference genome is very important, as it will be used as the backbone for the graph.  It is the only genome that is guaranteed to not have any cycles nor to ever be clipped, and therefore provides a coordinate system in the graph.  Any input genome can be used as a reference, provided it's consistently passed as the `--reference` option to all the commands.  It also must not have a "." in its genome name.  In practice, there are usually two possible references for the HPRC graphs: GRCh38 and CHM13.  To make a CHM13 graph, the process is identical to above, except the genome names should be renamed from `GRCh38 and CHM13_1_1.0` to `GRCh38.0 and CHM13_1_1` in the seqFile.  All `--reference GRCh38` flags need to be changed to `--reference CHM13_1_1`.

It is advisable to also pass `--vcfRerefence GRCh38` to `cactus-graphmap-join` to tell it to make a second VCF based on the GRCh38 reference. Likewise if creating a filtered graph with `cactus-graphmap-join --vgClipOpts "-d 9"`, you should use `--vgClipOpts "-d 9 -e GRCh38"` to not filter any nodes on the GRCh38 reference paths -- this will make it easier to use this graph (via surjection) to call variants on GRCh38 as well as CHM13.

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

**cactus-graphmap-join**

If the softmasked sequence has been kept unaligned with the various options described, it can be clipped out with the `--clipLength 1000000` option in `cactus-graphmap-join.  To clip based on BED files (ie those from `cactus-preprocess`) use the `--clipBed` option instead.

If the sequences were clipped instead of softmasked in `cactus-preprocess`, the HAL file will contain the chopped up sequences (with their fragment offsets described using suffixes of their names). Most downstream tools (CAT, Assembly Hubs) will not handle these well, so it's best to stitch the paths back together (leaving the clipped-out sequence unaligned).  This can be done via the `--unclipSeqFile` option.  


