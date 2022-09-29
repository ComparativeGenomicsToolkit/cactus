# Minigraph-Cactus Paper Methods

## HPRC Graphs

### Construction

The graphs were constructed as described [here](../pangenome.md#hprc-version-1.0-graphs).  They are available to download [here](https://github.com/human-pangenomics/hpp_pangenome_resources/).

The versions of the graph used to create Supplementary Figure 3 in order to show sequence removed when not using dna-brnn were created with the commands listed [here](../pangenome.md#hprc-graph).

### Graph Properties

#### Non-reference sequence

This was computed using [count-vg-hap-cov](https://github.com/ComparativeGenomicsToolkit/hal2vg/blob/f3d9a1838d1fb5582b6e1cd509792daee51fd2a9/count-vg-hap-cov.cpp) on the [GRCh38](https://s3-us-west-2.amazonaws.com/human-pangenomics/index.html?prefix=pangenomes/scratch/2021_08_11_minigraph_cactus/GRCh38-chrom-graphs/) and [CHM13](https://s3-us-west-2.amazonaws.com/human-pangenomics/index.html?prefix=pangenomes/scratch/2021_08_11_minigraph_cactus/CHM13-chrom-graphs/) chromosome vg graphs.


## D. Melanogaster Graphs

### Construction

The graphs and indexes were constructed using [these commands](../mc-pangenomes/16-fly-pg-2022-05-26-commands.md) with [this seqfile](../mc-pangenomes/16-fly-pg-2022-05-26-seqfile.txt) with Cactus commit [3f60d4f247c62d499e17202e059ff4d5d19cc71d](https://github.com/ComparativeGenomicsToolkit/cactus/commit/3f60d4f247c62d499e17202e059ff4d5d19cc71d).  The results are found at s3://vg-k8s/users/hickey/fly-pangenome/

The Progressive Cactus graph was created with Cactus commit [c31e4104b28ef8f55e2244177c593b9c2d7c39ac](https://github.com/ComparativeGenomicsToolkit/cactus/commit/c31e4104b28ef8f55e2244177c593b9c2d7c39ac) and default parameters

```
cactus ./js ./c31e4104b28ef8f55e2244177c593b9c2d7c39ac ./softmask.fly.hal --realTimeLogging --logFile softmask.fly.hal.log
hal2vg ./softmask.fly.hal --noAncesstors > fly-cactus.vg
```
and the output is available at s3://vg-k8s/users/hickey/fly-cactus/


### Short Read Sequencing Data

The reads were downloaded from NCBI as follows, using accessions listed in [fly-sra.tsv](./fly-sra.tsv).

```
for i in `cat fly-sra.tsv  | awk '{print $2}'` ; do fasterq-dump $i --split-3 ; bgzip ${i}_1.fastq --threads 8 ; bgzip ${i}_2.fastq --threads 8; done
```

### Mapping

Reads for accession SRR834526 (the same process was repeated for each sample) were mapped to the pangenome:

```
vg giraffe -m ./join-fly-pg-may26-d2/fly-pg-may26-d2.min -d ./join-fly-pg-may26-d2/fly-pg-may26-d2.dist -x ./join-fly-pg-may26-d2/fly-pg-may26-d2.xg -g ./join-fly-pg-may26-d2/fly-pg-may26-d2.gg -H ./join-fly-pg-may26-d2/fly-pg-may26-d2.gbwt -f  SRR834526_1.fastq.gz -f SRR834526_2.fastq.gz  -N $i -o gaf -t 32 > fly-pg-may26-d2_SRR834526.gaf
```

and mapped with BWA-MEM (v0.7.17)

```
bwa mem dm6.fa SRR834526_1.fastq.gz SRR834526_2.fastq.gz -t 16 > SRR834526_bwa.sam
```

### SV Genotyping

To genotype variants in the pangenome, vg pack computes the coverage of short reads in each node and edge of the pangenome, and vg call identifies the most likely paths supported by the short reads, effectively predicting genotypes at each snarl. The following commands were used:

```
vg pack -x  ./join-fly-pg-may26-d2/fly-pg-may26-d2.xg -a fly-pg-may26-d2_SRR834526.gaf.gz -Q 5 -o fly-pg-may26-d2_SRR834526.pack -t 16 
vg call ./join-fly-pg-may26-d2/fly-pg-may26-d2.xg -r./join-fly-pg-may26-d2/fly-pg-may26-d2.snarls -k  fly-pg-may26-d2_SRR834526.pack -a -A -t 16 -s SRR834526
```

For each sample, these variant calls were decomposed into canonical SVs using the same approach described above on the deconstructed VCF. The SV calls were then compared to the SVs in the pangenome using the sveval package [https://doi.org/10.1186/s13059-020-1941-7] which matches SVs based on their types, sizes and location. Here, two SVs were matched if: their regions had a reciprocal overlap of at least 90% for deletions and inversions; they were located at less than 100bp from each other, and their inserted sequences were at least 90% similar for insertions. The same approach was used to cluster the SVs alleles into the SV sites reported in the text and figures. The SV alleles were annotated with RepeatMasker (v4.0.9). We assigned a repeat class to a SV if more than 80% of the allelic sequence was annotated as such.

### Small Variants

We used `vg surject` to produce BAM files referenced on dm6 from the mappings to the pangenome, and FreeBayes v1.3.6 [https://doi.org/10.48550/arXiv.1207.3907] (in the absence of a high quality DeepVariant model) to call variants on these mappings and those from BWA-MEM. 

```
vg surject -x fly-pg-may26-d2.xg -t 16 -i -G fly-pg-may26-d2_SRR834526.gaf.gz -b -N $i | samtoolv view -O BAM > SRR834526_bwa.sam
```

and called with Freebayes (v1.3.6):

```
samtools sort fly-pg-may26-d2_SRR834526.bam --threads 2  > SRR834526.bam ; mv SRR834526.bam fly-pg-may26-d2_SRR834526.bam
freebayes -f dm6.fa  fly-pg-may26-d2_SRR834526.bam -=
```

```
samtools view -h SRR834526_bwa.sam -O BAM | samtools sort –threads 2 > SRR834526_bwa.bam
freebayes -f dm6.fa SRR834526_bwa.bam -=
```

Single-sample VCFs were merged with `bcftools merge`.

To compare the variant calls by both approaches, we used bcftools [https://doi.org/10.1093/gigascience/giab008] (v1.10.2) to normalize the VCFs (bcftools norm), and compare them (bcftools isec) to mark variant sites where both approaches call a variant, and sites where only one approach does. We compared the number of calls in each category, across samples, and for different minimum variant quality thresholds (QUAL field or genotype quality GQ field). 


