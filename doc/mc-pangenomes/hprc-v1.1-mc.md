# HPRC version 1.1 Minigraph-Cactus Release

This document describes the tools and commands used to make the version 1.1 Minigraph-Cactus HPRC graphs and indexes.

The main differences in 1.1 with respect to 1.0 are:

* Simpler, easier-to-reproduce pipeline (all steps to reproduce from scratch below)
* Output compatible with [current version of vg](https://github.com/vgteam/vg/releases/tag/v1.49.0), whose formats have changed considerably in many cases
* Cleaner base alignment via `minigraph`'s new base alignment options, as well as many bugfixes to `cactus`'s base aligner
* Cleaner large-scale alignment via better chaining and filtering.

For a list of all changes, consult the Release Notes for Cactus [version 2.2.1](https://github.com/ComparativeGenomicsToolkit/cactus/releases/tag/v2.2.1) through [version 2.6.4](https://github.com/ComparativeGenomicsToolkit/cactus/releases/tag/v2.6.4)

## Table of Contents

* [Graphs and Indexes](#graphs-and-indexes)
    * [Prepare the Input Assembly List](#prepare-the-input-assembly-list)
    * [Run cactus-pangenome to Make the Graphs and Indexes](#run-cactus-pangenome-to-make-the-graphs-and-indexes)
* [VCF Postprocessing](#vcf-postprocessing)
    * [VCFWave Decomposition](#vcfwave-decomposition)
    * [PanGenie Filtering](#pangenie-filtering)
    * [VCF Overview](#vcf-overview)
* [Excluded Regions](#excluded-regions)
* [MAF and TAF](#maf-and-taf)

## Graphs and Indexes

### Prepare the Input Assembly List

We download the assembly list from the HPRC website, and finesse it into a cactus seqfile

```
wget -q https://raw.githubusercontent.com/human-pangenomics/HPP_Year1_Assemblies/main/assembly_index/Year1_assemblies_v2_genbank.index
grep GRCh38 Year1_assemblies_v2_genbank.index | sed -e 's/_no_alt_analysis_set\t/\t/g' | awk '{print $1 "\t" $2}' > hprc-v1.1-mc.seqfile
tail -n +2 Year1_assemblies_v2_genbank.index | awk '{print $1 ".1\t" $2}' | grep -v CHM13 | grep -v GRCh38 >> hprc-v1.1-mc.seqfile
tail -n +2 Year1_assemblies_v2_genbank.index | awk '{print $1 ".2\t" $3}' | grep -v CHM13 | grep -v GRCh38 >> hprc-v1.1-mc.seqfile
sort -k1 hprc-v1.1-mc.seqfile > hprc-v1.1-mc.seqfile.sort ; mv hprc-v1.1-mc.seqfile.sort hprc-v1.1-mc.seqfile
sed hprc-v1.1-mc.seqfile -i -e 's%s3://human-pangenomics/working/%https://s3-us-west-2.amazonaws.com/human-pangenomics/working/%g'
```

For Year 1, we used CHM13 v1.1 with chrY from GRCh38 added (in future releases, we wil switch to CHM13 v2 which has chrY from HG002). Note that in the GRCh38-reference v1.0 graphs, we used a CHM13 without the Y, but we use the same for both graphs here to be consistent.
```
printf "CHM13\tftp://ftp.dfci.harvard.edu/pub/hli/minigraph/HPRC-f1g/CHM13v11Y.fa.gz\n" >> hprc-v1.1-mc.seqfile
```

We decided to leave three samples out for evaluation (if you want to make a graph with these samples included, just don't run this command)

```
grep -v 'HG002\|HG005\|NA19240' hprc-v1.1-mc.seqfile > t && mv t hprc-v1.1-mc.seqfile
```

Also, a misjoin in `HG02080#1#JAHEOW010000073.1` was manually corrected by using `samtools faidx` to break it into `HG02080#1#JAHEOW010000073.1_sub_0_7238466` and `HG02080#1#JAHEOW010000073.1_sub_7238466_12869124`. The `sub_X_Y` (0-based, open-ended like BED) coordinates are understood by the pipeline, and the offsets will be preserved in the GFA W-lines at the end. If we don't apply this change, then path names with ":"'s will end up in the HAL which will prevent it from working with assembly hubs.

```
wget -q $(grep HG02080\.1 hprc-v1.1-mc.seqfile | tail -1 | awk '{print $2}') -O HG02080.1.fa.gz
gzip -d HG02080.1.fa.gz
samtools faidx HG02080.1.fa
keep_contigs=$(awk '{print $1}' HG02080.1.fa.fai | grep -v JAHEOW010000073\.1)
samtools faidx HG02080.1.fa ${keep_contigs} > HG02080.1.fix.fa
samtools faidx HG02080.1.fa "HG02080#1#JAHEOW010000073.1:1-7238466" | sed -e 's/\([^:]*\):\([0-9]*\)-\([0-9]*\)/echo "\1_sub_$((\2-1))_\3"/e' >> HG02080.1.fix.fa
samtools faidx HG02080.1.fa "HG02080#1#JAHEOW010000073.1:7238467-12869124" | sed -e 's/\([^:]*\):\([0-9]*\)-\([0-9]*\)/echo "\1_sub_$((\2-1))_\3"/e' >> HG02080.1.fix.fa
bgzip HG02080.1.fix.fa --threads 8
```

Update the seqfile with the fixed assembly
```
grep -v HG02080\.1 hprc-v1.1-mc.seqfile > t && mv t hprc-v1.1-mc.seqfile
printf "HG02080.1\t./HG02080.1.fix.fa.gz\n" >> hprc-v1.1-mc.seqfile
```

### Run cactus-pangenome to Make the Graphs and Indexes

Make the GRCh38-based pangenome (run on SLURM head node) with [Cactus v2.6.4 binary release](https://github.com/ComparativeGenomicsToolkit/cactus/releases/tag/v2.6.4).  It took 31 hours wall time (details in logs, including memory use for each command)
```
cactus-pangenome ./js-grch38 ./hprc-v1.1-mc.seqfile --outName hprc-v1.1-mc-grch38 --outDir hprc-v1.1-mc-grch38 --reference GRCh38 CHM13 --filter 9 --giraffe clip filter --vcf  --viz --odgi --chrom-vg clip filter --chrom-og --gbz clip filter full --gfa clip full --vcf --logFile hprc-v1.1-mc-grch38.log --batchSystem slurm  --mgCores 64 --mapCores 16 --consCores 64 --indexCores 64 2> hprc-v1.1-mc-grch38.stderr &
```
And the CHM12-based pangenome (same command, just with reference order flipped and extra VCF specified), which took 22 hours wall time.
```
cactus-pangenome ./js-chm13 ./hprc-v1.1-mc.seqfile --outName hprc-v1.1-mc-chm13 --outDir hprc-v1.1-mc-chm13 --reference CHM13 GRCh38 --filter 9 --giraffe clip filter --vcf --vcfReference CHM13 GRCh38 --viz --odgi --chrom-vg clip filter --chrom-og --gbz clip filter full --gfa clip full --vcf --logFile hprc-v1.1-mc-chm13.log --batchSystem slurm  --mgCores 64 --mapCores 16 --consCores 64 --indexCores 64 2> hprc-v1.1-mc-chm13.stderr &
```

## VCF Postprocessing

### VCFWave Decomposition

This script reproduces the filtering used for the "decomposed" 1.0 VCFs, ex `hprc-v1.0-mc.grch38.vcfbub.a100k.wave.vcf.gz`.

It runs `vcfbub` to remove nested sites, keeping the highest-level version whose maximum allele length is under `100kb`, then runs `vcfwave` to re-align the alt alleles to the reference alleles, which is useful for cleaning up big muli-base SNPs, for example, by splitting them into smaller variants.  We use [vcflib commit b5287d1c2dbbfe53fc7af8f88ff4aa1208d5dd49](https://github.com/vcflib/vcflib/commits/b5287d1c2dbbfe53fc7af8f88ff4aa1208d5dd49) for `vcfwave` and the version of `vcfbub` from [Cactus v2.6.4 binary release](https://github.com/ComparativeGenomicsToolkit/cactus/releases/tag/v2.6.4).

This relies on the following script, `parallel-vcfwave.sh`

```
#!/bin/bash
INFILE=$1
OUTFILE=$2
for chr in `bcftools view -h $INFILE | perl -ne 'if (/^##contig=<ID=([^,]+)/) { print "$1\n" }'`; do echo $chr; done | parallel -j 16 "bcftools view $INFILE -r {} | bgzip > $INFILE-{}.input.vcf.gz ; tabix -fp vcf $INFILE-{}.input.vcf.gz ; vcfbub --input $INFILE-{}.input.vcf.gz -l 0 -a 100000 | vcfwave -t 2 -I 1000 | bgzip --threads 2 > ${INFILE}-{}.vcf.gz; rm -f $INFILE-{}.input.vcf.gz $INFILE-{}.input.vcf.gz.tbi"
bcftools concat ${INFILE}-*.vcf.gz | bcftools sort | bgzip --threads 16 > ${OUTFILE}
tabix -fp vcf ${OUTFILE}
rm -f ${INFILE}-*.vcf.gz
```

```
./parallel-vcfwave.sh hprc-v1.1-mc-grch38.raw.vcf.gz hprc-v1.1-mc-grch38.vcfbub.a100k.wave.vcf.gz
./parallel-vcfwave.sh hprc-v1.1-mc-chm13.raw.vcf.gz hprc-v1.1-mc-chm13.vcfbub.a100k.wave.vcf.gz
./parallel-vcfwave.sh hprc-v1.1-mc-chm13.GRCh38.raw.vcf.gz hprc-v1.1-mc-chm13.GRCh38.vcfbub.a100k.wave.vcf.gz
```

### VCF Overview

Three distinct VCFs were produced from the graphs.

* `hprc-v1.1-mc-grch38.raw.vcf.gz` : GRCh38-based VCF from the GRCh38-based graph, including nested variants
* `hprc-v1.1-mc-chm13.raw.vcf.gz` : CHM13-based VCF from the CHM13-based graph, including nested variants
* `hprc-v1.1-mc-chm13.GRCh38.raw.vcf.gz` : GRCh38-based VCF from the CHM13-based graph, including nested variants

All three VCFs are filtered independently, in the same way, so below we will only discuss `hprc-v1.1-mc-grch38.raw.vcf.gz`.

`minigraph-cactus` also outputs by default:

* `hprc-v1.1-mc-grch38.vcf.gz`. Nested sites, and those with reference allele length > 100kb are removed with `vcfbub --max-ref-length 100000 --max-level 0`.  Note that it can still contain alt alleles > 100kb. The sites in this VCF are consistent with the graph topoloy, with the alleles described as paths in the (GFA) graph via their `AT` tags. 

The more aggressive filtering procedure, described above, was also run similar to the 1.0 results, where each graph was derived from the previous one:

1) `hprc-v1.1-mc-grch38.raw.vcf.gz` : Raw output, as described above.
2) `hprc-v1.1-mc-grch38.vcfbub.a100k.wave.vcf.gz` : Variants decomposed with `vcfbub -l 0 -a 100000` and re-aligned with `vcfwave -I 1000`.

**Important** Because of the `vcfwave` realignment, the alleles in `hprc-v1.1-mc-grch38.vcfbub.a100k.wave.vcf.gz` are no longer necessarily consistent with the topology of the graph.


## Excluded Regions

Using `samtools faidx` on each contig in the seqfile, we can obtain a [BED file of all input contigs](hprc-v1.1-mc-input-contigs.bed).  The path fragments in the graphs, as obtained form the GFA W-lines, can then be subtracted from this file in order to see what's missing from the graphs.

`extract-missing.sh`
```
#!/bin/bash
INGFA=$1
zcat $INGFA | grep ^W | awk '{print $2 "#" $3 "#" $4 "\t" $5 "\t" $6}' | bedtools sort > ${INGFA}.contigs.bed
bedtools subtract -a hprc-v1.1-mc-input-contigs.bed -b ${INGFA}.contigs.bed 
```

```
./extract-missing.sh hprc-v1.1-mc-grch38.full.gfa.gz > hprc-v1.1-mc-grch38.full.cut.bed &
./extract-missing.sh hprc-v1.1-mc-grch38.gfa.gz > hprc-v1.1-mc-grch38.cut.bed &

./extract-missing.sh hprc-v1.1-mc-chm13.full.gfa.gz > hprc-v1.1-mc-chm13.full.cut.bed &
./extract-missing.sh hprc-v1.1-mc-chm13.gfa.gz > hprc-v1.1-mc-chm13.cut.bed &

```

The `.full.cut.bed` files will show contigs that could not be assigned to a chromosome, and therefore are cut from the full graph. The `.cut.bed` files will also include contig fragments removed when clipping out unaligned regions >10kb.  

## MAF and TAF

[Multiple Alignment Format (MAF)](https://genome.ucsc.edu/FAQ/FAQformat.html#format5) is a text-based format that encodes alignment blocks. [Transposed Alignment Format (TAF)](https://github.com/ComparativeGenomicsToolkit/taffy#the-taf-file-format) is also a test-based multiple alignment format but whereas MAF stores one row per line, TAF stores one column per line.  The two formats encode the same information, but TAF tends to be much smaller as it can use run-length encoding to compress low-information columns, and is more clever about not repeating sequence names every line.

When converting from HAL to MAF (or TAF), the alignment is projected onto a chosen reference (we use GRCh38 and CHM13 for the GRCh38-based and CHM13-based graphs, respectively). Alignments that cannot be represented within the context of the chosen reference are not included (so the HAL->MAF/TAF conversion is lossy).

Both MAF and TAF files are indexed for random access via reference coordinates.  The indexes are `.tai` files, and work similarly to `samtools` FASTA indexes.  They can be queried with `taffy view`.  The `taffy` binary is included with Cactus and its source code and documentation can be found [here](https://github.com/ComparativeGenomicsToolkit/taffy). `taffy` works with both MAF and TAF.

The HAL alignments were converted to MAF as follows:

First, the genomes were renamed so that "." characters in the sample names were replaced by "#" characters (since "."s take on special meaning in some MAF tools including `mafDuplicateFilter` below).  This is something that needs to be smoothed over in future versions of the toolset.  See below for the contents of `hprc-rename.tsv`.

```
halRenameGenome hprc-v1.1-mc-grch38.full.hal hprc-rename.tsv
halRenameGenome hprc-v1.1-mc-chm13.full.hal hprc-rename.tsv
```

Next, the HAL files are converted to MAF (Cactus version, which includes all binaries referred to below, used: v2.6.7).


```
cactus-hal2maf ./js ./hprc-v1.1-mc-grch38.full.hal ./hprc-v1.1-mc-grch38.full.maf.gz --noAncestors --refGenome GRCh38 --filterGapCausingDupes --chunkSize 100000 --batchCores 96 --batchCount 10 --noAncestors --batchParallelTaf 32 --batchSystem slurm --logFile ./hprc-v1.1-mc-grch38.full.maf.gz.log

cactus-hal2maf ./js ./hprc-v1.1-mc-chm13.full.hal ./hprc-v1.1-mc-chm13.full.maf.gz --noAncestors --refGenome CHM13 --filterGapCausingDupes --chunkSize 100000 --batchCores 96 --batchCount 10 --noAncestors --batchParallelTaf 32 --batchSystem slurm --logFile ./hprc-v1.1-mc-chm13.full.maf.gz.log

```

Duplications (blocks with multiple rows from the same haplotype) are filtered out, as some tools do not expect them:

```
zcat  ./hprc-v1.1-mc-grch38.full.maf.gz | mafDuplicateFilter -k -m - | bgzip > hprc-v1.1-mc-grch38.full.single-copy.maf.gz
zcat  ./hprc-v1.1-mc-chm13.full.maf.gz | mafDuplicateFilter -k -m - | bgzip > hprc-v1.1-mc-chm13.full.single-copy.maf.gz
```

Index the MAF files
```
taffy index -i hprc-v1.1-mc-grch38.full.maf.gz
taffy index -i hprc-v1.1-mc-chm13.full.maf.gz
taffy index -i hprc-v1.1-mc-grch38.full.single-copy.maf.gz
taffy index -i hprc-v1.1-mc-chm13.full.single-copy.maf.gz
```

Create and index the TAF versions (they are way smaller)
```
taffy view -i hprc-v1.1-mc-grch38.full.maf.gz -uc > hprc-v1.1-mc-grch38.full.taf.gz
taffy view -i hprc-v1.1-mc-chm13.full.maf.gz -uc > hprc-v1.1-mc-chm13.full.taf.gz
taffy view -i hprc-v1.1-mc-grch38.full.single-copy.maf.gz -uc > hprc-v1.1-mc-grch38.full.single-copy.taf.gz
taffy view -i hprc-v1.1-mc-chm13.full.single-copy.maf.gz -uc > hprc-v1.1-mc-chm13.full.single-copy.taf.gz

taffy index -i hprc-v1.1-mc-grch38.full.taf.gz
taffy index -i hprc-v1.1-mc-chm13.full.taf.gz
taffy index -i hprc-v1.1-mc-grch38.full.single-copy.taf.gz
taffy index -i hprc-v1.1-mc-chm13.full.single-copy.taf.gz

```



The contenst of `hprc-rename.tsv`

```
NA20129.2       NA20129#2
NA21309.2       NA21309#2
HG02572.1       HG02572#1
NA18906.2       NA18906#2
NA18906.1       NA18906#1
HG03492.2       HG03492#2
HG03492.1       HG03492#1
HG03486.2       HG03486#2
HG03486.1       HG03486#1
HG02148.2       HG02148#2
HG01361.2       HG01361#2
HG01928.2       HG01928#2
HG02257.1       HG02257#1
HG01891.1       HG01891#1
HG02486.1       HG02486#1
HG02257.2       HG02257#2
HG01358.2       HG01358#2
NA21309.1       NA21309#1
HG01106.1       HG01106#1
HG01358.1       HG01358#1
HG01258.2       HG01258#2
HG01258.1       HG01258#1
HG00741.1       HG00741#1
HG01952.2       HG01952#2
HG00673.2       HG00673#2
HG01175.1       HG01175#1
HG01071.2       HG01071#2
HG00621.1       HG00621#1
HG02148.1       HG02148#1
HG01123.2       HG01123#2
HG00621.2       HG00621#2
HG00438.1       HG00438#1
HG03453.2       HG03453#2
HG02723.2       HG02723#2
HG00438.2       HG00438#2
HG00673.1       HG00673#1
HG01175.2       HG01175#2
HG02486.2       HG02486#2
HG00735.1       HG00735#1
HG02717.1       HG02717#1
HG01071.1       HG01071#1
HG02145.1       HG02145#1
HG00741.2       HG00741#2
HG01361.1       HG01361#1
HG02622.2       HG02622#2
HG01978.2       HG01978#2
HG01106.2       HG01106#2
HG02055.2       HG02055#2
HG02109.1       HG02109#1
HG01928.1       HG01928#1
HG01123.1       HG01123#1
HG03579.1       HG03579#1
HG02559.1       HG02559#1
HG02559.2       HG02559#2
HG01952.1       HG01952#1
HG02572.2       HG02572#2
HG01978.1       HG01978#1
HG03579.2       HG03579#2
HG02622.1       HG02622#1
HG02630.1       HG02630#1
HG02630.2       HG02630#2
HG02080.2       HG02080#2
HG02818.2       HG02818#2
HG02717.2       HG02717#2
HG02886.1       HG02886#1
HG02886.2       HG02886#2
HG02145.2       HG02145#2
HG01891.2       HG01891#2
HG03540.1       HG03540#1
HG03516.1       HG03516#1
HG03516.2       HG03516#2
HG03540.2       HG03540#2
HG01243.1       HG01243#1
HG00733.1       HG00733#1
HG02818.1       HG02818#1
HG00735.2       HG00735#2
HG00733.2       HG00733#2
HG01109.1       HG01109#1
HG01109.2       HG01109#2
HG03453.1       HG03453#1
HG02723.1       HG02723#1
HG02055.1       HG02055#1
HG01243.2       HG01243#2
HG02080.1       HG02080#1
HG02109.2       HG02109#2
NA20129.1       NA20129#1
HG03098.1       HG03098#1
HG03098.2       HG03098#2
```