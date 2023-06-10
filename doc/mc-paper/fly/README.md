# D. Melanogaster Pangenome Methods

## Construction

The graphs and indexes were constructed using [these commands](../../mc-pangenomes/16-fly-pg-2022-05-26-commands.md) with [this seqfile](../../mc-pangenomes/16-fly-pg-2022-05-26-seqfile.txt) with Cactus commit [3f60d4f247c62d499e17202e059ff4d5d19cc71d](https://github.com/ComparativeGenomicsToolkit/cactus/commit/3f60d4f247c62d499e17202e059ff4d5d19cc71d).  The results (and an archive of the input fastas) are found [here](https://s3-us-west-2.amazonaws.com/human-pangenomics/index.html?prefix=publications/mc_2022/mc_pangenomes/16-fruitfly-mc-2022-05-26/).

The Progressive Cactus graph was created with Cactus commit [c31e4104b28ef8f55e2244177c593b9c2d7c39ac](https://github.com/ComparativeGenomicsToolkit/cactus/commit/c31e4104b28ef8f55e2244177c593b9c2d7c39ac) and default parameters

```
cactus ./js 16-fly-pg-2022-05-26-seqfile.txt ./softmask.fly.hal --realTimeLogging --logFile softmask.fly.hal.log
hal2vg ./softmask.fly.hal --noAncesstors > fly-cactus.vg
```

and the output is available [here](https://s3-us-west-2.amazonaws.com/human-pangenomics/index.html?prefix=publications/mc_2022/fruitfly/fly-cactus/).


## Short Read Sequencing Data

The reads were downloaded from NCBI as follows, using accessions listed in [fly-sra.tsv](../fly-sra.tsv) and `faster-dump` from [sra toolkit v3.0.2](https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/3.0.2/sratoolkit.3.0.2-ubuntu64.tar.gz).

```
for i in `cat fly-sra.tsv  | awk '{print $2}'` ; do fasterq-dump $i --split-3 ; bgzip ${i}_1.fastq --threads 8 ; bgzip ${i}_2.fastq --threads 8; done
```

## Mapping

Reads for accession SRR834526 (the same process was repeated for each sample) were mapped to the pangenome:

```
vg giraffe -m 16-fruitfly-mc-2022-05-26-d2.min -d 16-fruitfly-mc-2022-05-26-d2.dist -x 16-fruitfly-mc-2022-05-26-d2.xg -g 16-fruitfly-mc-2022-05-26-d2.gg -H 16-fruitfly-mc-2022-05-26-d2.gbwt -f  SRR834526_1.fastq.gz -f SRR834526_2.fastq.gz  -N $i -o gaf -t 32 > 16-fruitfly-mc-2022-05-26-d2_SRR834526.gaf
```

and mapped with BWA-MEM (v0.7.17)

```
bwa mem dm6.fa SRR834526_1.fastq.gz SRR834526_2.fastq.gz -t 16 > SRR834526_bwa.sam
```

## SV Genotyping

To genotype variants in the pangenome, vg pack computes the coverage of short reads in each node and edge of the pangenome, and vg call identifies the most likely paths supported by the short reads, effectively predicting genotypes at each snarl. The following commands were used:

```
vg pack -x  16-fruitfly-mc-2022-05-26-d2.xg -a 16-fruitfly-mc-2022-05-26-d2_SRR834526.gaf.gz -Q 5 -o 16-fruitfly-mc-2022-05-26-d2_SRR834526.pack -t 16 
vg call 16-fruitfly-mc-2022-05-26-d2.xg -r16-fruitfly-mc-2022-05-26-d2.snarls -k  16-fruitfly-mc-2022-05-26-d2_SRR834526.pack -a -A -t 16 -s SRR834526 > 16-fruitfly-mc-2022-05-26-d2_SRR834526.vcf
```

For each sample, these variant calls were decomposed into canonical SVs using the same approach described above on the deconstructed VCF. The SV calls were then compared to the SVs in the pangenome using the sveval package [https://doi.org/10.1186/s13059-020-1941-7] which matches SVs based on their types, sizes and location. Here, two SVs were matched if: their regions had a reciprocal overlap of at least 90% for deletions and inversions; they were located at less than 100bp from each other, and their inserted sequences were at least 90% similar for insertions. The same approach was used to cluster the SVs alleles into the SV sites reported in the text and figures. The SV alleles were annotated with RepeatMasker (v4.0.9). We assigned a repeat class to a SV if more than 80% of the allelic sequence was annotated as such.

## Small Variants

We used `vg surject` to produce BAM files referenced on dm6 from the mappings to the pangenome, and FreeBayes v1.3.6 [https://doi.org/10.48550/arXiv.1207.3907] (in the absence of a high quality DeepVariant model) to call variants on these mappings and those from BWA-MEM. 

```
vg surject -x 16-fruitfly-mc-2022-05-26-d2.xg -t 16 -i -G 16-fruitfly-mc-2022-05-26-d2_SRR834526.gaf.gz -b -N $i | samtoolv view -O BAM > SRR834526_bwa.sam
```

and called with Freebayes (v1.3.6):

```
samtools sort 16-fruitfly-mc-2022-05-26-d2_SRR834526.bam --threads 2  > SRR834526.bam ; mv SRR834526.bam 16-fruitfly-mc-2022-05-26-d2_SRR834526.bam
freebayes -f dm6.fa  16-fruitfly-mc-2022-05-26-d2_SRR834526.bam -=
```

```
samtools view -h SRR834526_bwa.sam -O BAM | samtools sort –threads 2 > SRR834526_bwa.bam
freebayes -f dm6.fa SRR834526_bwa.bam -=
```

Single-sample VCFs were merged with `bcftools merge`.

To compare the variant calls by both approaches, we used bcftools [https://doi.org/10.1093/gigascience/giab008] (v1.10.2) to normalize the VCFs (bcftools norm), and compare them (bcftools isec) to mark variant sites where both approaches call a variant, and sites where only one approach does. We compared the number of calls in each category, across samples, and for different minimum variant quality thresholds (QUAL field or genotype quality GQ field). 



## Exploration of the mapping and variant calls

### Snakemake pipeline, scripts, and dependencies

[Snakemake](https://snakemake.github.io/) is used for multiple analysis, to help run the same analysis across many samples (e.g. compute the mapping statistics across the 100 samples).
For this analysis the raw data includes: 
- the pangenome indexes (full graph and frequency-filtered graph *-d2*)
- the deconstructed VCF
- the alignments of the 100 samples to the pangenome (GAMs) and to dm6 (gzipped SAM files)
- the variants genotyped from the alignment to the pangenome for the 100 samples
- the calls from FreeBayes for the 100 samples for both the pangenomic and linear reference approaches.

To reproduce the analysis or generate the commands used, the data could be downloaded/mimicked locally using the following commands:

```sh
mkdir -p pangenomes reads calls vg_calls
## pangenomes and deconstructed VCF
wget -O pangenomes/16-fruitfly-mc-2022-05-26.xg https://s3-us-west-2.amazonaws.com/human-pangenomics/publications/mc_2022/mc_pangenomes/16-fruitfly-mc-2022-05-26/16-fruitfly-mc-2022-05-26.xg
wget -O pangenomes/16-fruitfly-mc-2022-05-26.vcf.gz https://s3-us-west-2.amazonaws.com/human-pangenomics/publications/mc_2022/mc_pangenomes/16-fruitfly-mc-2022-05-26/16-fruitfly-mc-2022-05-26.vcf.gz
wget -O pangenomes/16-fruitfly-mc-2022-05-26-d2.xg https://s3-us-west-2.amazonaws.com/human-pangenomics/publications/mc_2022/mc_pangenomes/16-fruitfly-mc-2022-05-26/16-fruitfly-mc-2022-05-26-d2.xg
## mapping and calls across 100 samples
for samp in SRR833244 SRR834507 SRR834508 SRR834509 SRR834510 SRR834511 SRR834512 SRR834514 SRR834516 SRR834517 SRR834519 SRR834520 SRR834521 SRR834522 SRR834523 SRR834526 SRR834527 SRR834531 SRR834537 SRR834539 SRR834541 SRR834543 SRR834545 SRR834546 SRR834547 SRR834551 SRR834552 SRR834553 SRR834554 SRR835023 SRR835024 SRR835025 SRR835026 SRR835027 SRR835028 SRR835029 SRR835030 SRR835031 SRR835033 SRR835034 SRR835035 SRR835036 SRR835037 SRR835038 SRR835039 SRR835040 SRR835041 SRR835042 SRR835043 SRR835044 SRR835045 SRR835046 SRR835047 SRR835048 SRR835050 SRR835051 SRR835052 SRR835054 SRR835055 SRR835058 SRR835059 SRR835060 SRR835061 SRR835062 SRR835063 SRR835067 SRR835069 SRR835071 SRR835072 SRR835073 SRR835074 SRR835075 SRR835077 SRR835086 SRR835087 SRR835088 SRR835089 SRR835091 SRR835092 SRR835095 SRR835096 SRR835097 SRR835098 SRR932121 SRR933563 SRR933569 SRR933573 SRR933577 SRR933580 SRR933581 SRR933585 SRR933586 SRR933587 SRR933589 SRR933591 SRR933592 SRR933593 SRR933594 SRR933599 SRR933601
do
    wget -O reads/16-fruitfly-mc-2022-05-26-d2.${samp}.giraffe.gaf.gz https://s3-us-west-2.amazonaws.com/human-pangenomics/publications/mc_2022/fruitfly/fly-mappings/fly-pg-may26-d2_${samp}.gaf.gz
    wget -O reads/${samp}_bwa.sam.gz https://s3-us-west-2.amazonaws.com/human-pangenomics/publications/mc_2022/fruitfly/fly-bwa/${samp}_bwa.sam.gz
    ## we shared the merged VCFs of these so they are not available because redundant
    ## "touch"ing them allows us to see how they were processed by the snakemake workflow
    touch calls/${samp}_bwa.fb.vcf.gz 
    touch calls/${samp}_bwa.fb.vcf.gz.tbi
    touch calls/16-fruitfly-mc-2022-05-26-d2_${samp}.fb.vcf.gz
    touch calls/16-fruitfly-mc-2022-05-26-d2_${samp}.fb.vcf.gz.tbi
    touch vg_calls/16-fruitfly-mc-2022-05-26-d2_${samp}.vcf.gz
done
```

We have also shared the log files enumerating all the commands ran by Snakemake (see `snakemake_log.*` files).

Scripts used by this workflow were deposited in the [`scripts`](scripts) folder.

Finally, we ran this analysis using a docker container with the necessary tools. 
It can be built using [this Dockerfile](Dockerfile).
In particular, the tools that were used in the analysis below include:

- [Snakemake](snakemake.readthedocs.io/)
- [bcftools](https://samtools.github.io/bcftools/bcftools.html)
- [vg](https://github.com/vgteam/vg)
- [jq](https://stedolan.github.io/jq/)
- [Python 3](https://www.python.org/) with the [cyvcf2](http://brentp.github.io/cyvcf2/) and [libbdsg](https://github.com/vgteam/libbdsg) (python binding for vg)
- [R](https://www.r-project.org/) with the [GenomicRanges](https://bioconductor.org/packages/release/bioc/html/GenomicRanges.html), [sveval](https://github.com/jmonlong/sveval), [dplyr](https://dplyr.tidyverse.org/), and [ggplot2](https://ggplot2.tidyverse.org/) packages

### Mapping statistics

To compute the mapping statistics from the GAMs (from aligning reads with Giraffe to the pangenome) and the gzipped SAMs (from aligning reads with BWA-MEM to dm6):

```sh
snakemake --cores 8 --config mapper=giraffe graph=16-fruitfly-mc-2022-05-26-d2 -p mapping_stats
snakemake --cores 8 --config mapper=bwa graph=dm6 -p mapping_stats
```

For detailed commands, see the log files of these two workflows: [snakemake_log.mapping_stats_mc_giraffe.log](snakemake_log.mapping_stats_mc_giraffe.log) and [snakemake_log.mapping_stats_dm6_bwa.log](snakemake_log.mapping_stats_dm6_bwa.log)

This will create a tsv for each sample and mapper in the forms: `<GENOME>_mappings/<GENOME>.<SAMPLE>.<MAPPER>.mapstats.txt` where:
- `<GENOME>` is either `dm6` or `16-fruitfly-mc-2022-05-26-d2`
- `<MAPPER>` is either `bwa` or `giraffe`

Those files, localized in the two `*_mappings` folders, are then read by [`mapstats-analysis.R`](mapstats-analysis.R).
This script will compute stats and makes the figures shown in the manuscript.
Those files were archived at [`https://s3-us-west-2.amazonaws.com/human-pangenomics/publications/mc_2022/fruitfly/fly-maptstats.tar.gz`](https://s3-us-west-2.amazonaws.com/human-pangenomics/publications/mc_2022/fruitfly/fly-maptstats.tar.gz).

### FreeBayes calls

To compare the small variant calls using bcftools:

```sh
snakemake --cores 8 --config graph=16-fruitfly-mc-2022-05-26-d2 -p compare_freebayes
```

For detailed commands, see the log file: [snakemake_log.compare_freebayes_calls.log](snakemake_log.compare_freebayes_calls.log)

This will create one `<SAMPLE>.isec.stats.tsv.gz` for each sample.
Assuming all those files are located in a `freebayes_stats` folder, the [`compare-freebayes-calls.R`](compare-freebayes-calls.R) script will compute stats and make the graphs shown in the manuscript.

---

To compare the proportion of heterozygous variants, the calls were first merged with:

```sh
snakemake --cores 8 -p results/bwa.dm6.vcf.gz
snakemake --cores 8 -p results/surject.16-fruitfly-mc-2022-05-26-d2.vcf.gz
```

For detailed commands, see the log files of these two workflows: [snakemake_log.merge_freebayes_mc_giraffe_calls.log](snakemake_log.merge_freebayes_mc_giraffe_calls.log) and [snakemake_log.merge_freebayes_dm6_bwa_calls.log](snakemake_log.merge_freebayes_dm6_bwa_calls.log).

The merged VCFs are made available at [`https://s3-us-west-2.amazonaws.com/human-pangenomics/publications/mc_2022/fruitfly/bwa.dm6.vcf.gz`](https://s3-us-west-2.amazonaws.com/human-pangenomics/publications/mc_2022/fruitfly/bwa.dm6.vcf.gz) and [`https://s3-us-west-2.amazonaws.com/human-pangenomics/publications/mc_2022/fruitfly/surject.16-fruitfly-mc-2022-05-26-d2.vcf.gz`](https://s3-us-west-2.amazonaws.com/human-pangenomics/publications/mc_2022/fruitfly/surject.16-fruitfly-mc-2022-05-26-d2.vcf.gz).

The [`compare-freebayes-calls-hetrate.R`](compare-freebayes-calls-hetrate.R) script then compute the heterozygous rates and makes the graph shown in the manuscript.
This info is also compiled in `freebayes-het-rate.tsv` which will be read by this script if present.
This file is available at [`https://s3-us-west-2.amazonaws.com/human-pangenomics/publications/mc_2022/fruitfly/freebayes-het-rate.tsv`](https://s3-us-west-2.amazonaws.com/human-pangenomics/publications/mc_2022/fruitfly/freebayes-het-rate.tsv)

### SV analysis

#### Decomposing variants

To decompose the deconstructed VCF from the full pangenome:

```sh
snakemake --cores 8 -p decompose_pangenome
```

For detailed commands, see the log file: [snakemake_log.decompose_pangenome.log](snakemake_log.decompose_pangenome.log).

This will create one file for each sample/assembly in the pangenome: `svs/<SAMPLE>.<GRAPH>.decomposed.svs.rds`

---

To decompose the calls for each of the 100 samples that were genotyped with the pangenome:

```sh
snakemake --cores 8 --config graph=16-fruitfly-mc-2022-05-26-d2 -p decompose_calls
```

For detailed commands, see the log file: [snakemake_log.decompose_vg_calls.log](snakemake_log.decompose_vg_calls.log).

The variants will be merged into one file called `16-fruitfly-mc-2022-05-26-d2.100samples.decomposed.svs.rds`

#### Post-decomposition/calling analysis

1. Cluster SVs into SV sites (for both SVs in the pangenome and calls across 100 samples): [`cluster-svs.R`](cluster-svs.R)
    - This will just add new columns to the calls to specify the assigned SV site.
    - Input: 
        - `svs/<SAMPLE>.<GRAPH>.decomposed.svs.rds` for each sample in the pangenome
        - `16-fruitfly-mc-2022-05-26-d2.100samples.decomposed.svs.rds`
    - Output:
        - `16-fruitfly-mc-2022-05-26.svs.site.rol90.insd100.rds`
        - `16-fruitfly-mc-2022-05-26-d2.100samples.decomposed.svs.site.rol90.insd100.rds`
    - The output files were deposited in:
        - [`https://s3-us-west-2.amazonaws.com/human-pangenomics/publications/mc_2022/fruitfly/16-fruitfly-mc-2022-05-26.svs.site.rol90.insd100.rds`](https://s3-us-west-2.amazonaws.com/human-pangenomics/publications/mc_2022/fruitfly/16-fruitfly-mc-2022-05-26.svs.site.rol90.insd100.rds)
        - [`https://s3-us-west-2.amazonaws.com/human-pangenomics/publications/mc_2022/fruitfly/16-fruitfly-mc-2022-05-26-d2.100samples.decomposed.svs.site.rol90.insd100.rds`](https://s3-us-west-2.amazonaws.com/human-pangenomics/publications/mc_2022/fruitfly/16-fruitfly-mc-2022-05-26-d2.100samples.decomposed.svs.site.rol90.insd100.rds)
2. Annotate the alleles with RepeatMasker: [`annotate-repeats.R`](annotate-repeats.R)
    - Annotates the major allele of each SV site
    - Uses the following docker container with RepeatMasker: `jmonlong/repeatmasker:release-4.0.9-p2`
    - Input: `16-fruitfly-mc-2022-05-26-d2.100samples.decomposed.svs.site.rol90.insd100.rds`
    - Output: `16-fruitfly-mc-2022-05-26-d2.100samples.decomposed.svs.site.rol90.insd100.rmsk.rds`
    - The output files were deposited at [`https://s3-us-west-2.amazonaws.com/human-pangenomics/publications/mc_2022/fruitfly/16-fruitfly-mc-2022-05-26-d2.100samples.decomposed.svs.site.rol90.insd100.rmsk.rds`](https://s3-us-west-2.amazonaws.com/human-pangenomics/publications/mc_2022/fruitfly/16-fruitfly-mc-2022-05-26-d2.100samples.decomposed.svs.site.rol90.insd100.rmsk.rds)
1. Make figures: [`svs-analysis.R`](svs-analysis.R).
    - Input:
        - `16-fruitfly-mc-2022-05-26.svs.site.rol90.insd100.rds`
        - `16-fruitfly-mc-2022-05-26-d2.100samples.decomposed.svs.site.rol90.insd100.rds`
        - `16-fruitfly-mc-2022-05-26-d2.100samples.decomposed.svs.site.rol90.insd100.rmsk.rds`

