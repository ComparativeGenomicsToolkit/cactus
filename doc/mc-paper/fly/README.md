# Analysis of the Fly pangenome

## Snakemake pipeline and scripts

[Snakemake](https://snakemake.github.io/) is used for multiple analysis, to help run the same analysis across many samples (e.g. compute the mapping statistics across the 100 samples).
Our pipeline, implemented in the [`Snakefile`](Snakefile), is set up to work with our data located in our private S3 buckets.
Hence, to reproduce the analysis, one would have to replace the parts specifying those S3 paths (`S3.remote(...)`), with local paths where the raw data was downloaded.
For this analysis the raw data includes: 
- the pangenome indexes (full graph and frequency-filtered graph *-d2*)
- the deconstructed VCF
- the alignments of the 100 samples to the pangenome (GAMs) and to dm6 (gzipped SAM files)
- the variants genotyped from the alignment to the pangenome for the 100 samples
- the calls from FreeBayes for the 100 samples for both the pangenomic and linear reference approaches.

The instructions to produce this input data can be found at ??.

The pipeline also uses custom scripts (`S3.remote(SROOT + 'scripts/...')` in the Snakefile). 
Those scripts were deposited in the [`scripts`](scripts) folder.

## Mapping statistics

To compute the mapping statistics from the GAMs (from aligning reads with Giraffe to the pangenome) and the gzipped SAMs (from aligning reads with BWA-MEM to dm6):

```sh
snakemake --cores 8 --config mapper=giraffe graph=fly-pg-may26-d2 -p mapping_stats
snakemake --cores 8 --config mapper=bwa graph=dm6 -p mapping_stats
```

This will create a tsv for each sample and mapper in the forms: `<GENOME>_mappings/<GENOME>.<SAMPLE>.<MAPPER>.mapstats.txt` where:
- *<GENOME>* is either `dm6` or `fly-pg-may26-d2`
- *<MAPPER>* is either `bwa` or `giraffe`

Those files, localized in the two `*_mappings` folders, are then read by [`mapstats-analysis.R`](mapstats-analysis.R).
This script will compute stats and makes the figures shown in the manuscript.
Those files were archived at [`http://public.gi.ucsc.edu/~jmonlong/hprc/minigraph-cactus-manuscript/fly/fly-maptstats.tar.gz`](http://public.gi.ucsc.edu/~jmonlong/hprc/minigraph-cactus-manuscript/fly/fly-maptstats.tar.gz).

## FreeBayes calls

To compare the small variant calls using bcftools:

```sh
snakemake --cores 8 --config graph=fly-pg-may26-d2 -p compare_freebayes
```

This will create one `<SAMPLE>.isec.stats.tsv.gz` for each sample.
Assuming all those files are in located in a `freebayes_stats` folder, the [`compare-freebayes-calls.R`](compare-freebayes-calls.R) script will compute stats and make the graphs shown in the manuscript.

---

To compare the proportion of heterozygous variants, the calls were first merged with:

```sh
snakemake --cores 8 -p bwa.dm6.vcf.gz
snakemake --cores 8 -p surject.fly-pg-may26-d2.vcf.gz
```

The merged VCFs are made available at [`http://public.gi.ucsc.edu/~jmonlong/hprc/minigraph-cactus-manuscript/fly/bwa.dm6.vcf.gz`](http://public.gi.ucsc.edu/~jmonlong/hprc/minigraph-cactus-manuscript/fly/bwa.dm6.vcf.gz) and [`http://public.gi.ucsc.edu/~jmonlong/hprc/minigraph-cactus-manuscript/fly/surject.fly-pg-may26-d2.vcf.gz`](http://public.gi.ucsc.edu/~jmonlong/hprc/minigraph-cactus-manuscript/fly/surject.fly-pg-may26-d2.vcf.gz).

The [`compare-freebayes-calls-hetrate.R`](compare-freebayes-calls-hetrate.R) script then compute the heterozygous rates and makes the graph shown in the manuscript.
This info is also compiled in `freebayes-het-rate.tsv` which will be read by this script if present.
This file is available at [`http://public.gi.ucsc.edu/~jmonlong/hprc/minigraph-cactus-manuscript/fly/freebayes-het-rate.tsv`](http://public.gi.ucsc.edu/~jmonlong/hprc/minigraph-cactus-manuscript/fly/freebayes-het-rate.tsv)

## SV analysis

### Decomposing variants

To decompose the deconstructed VCF from the full pangenome:

```sh
snakemake --cores 8 -p decompose_pangenome
```

This will create one file for each sample/assembly in the pangenome: `<SAMPLE>.decomposed.svs.rds`

---

To decompose the calls for each of the 100 samples that were genotyped with the pangenome:

```sh
snakemake --cores 8 --config graph=fly-pg-may26-d2 -p decompose_calls
```

The variants will be merged into one file called `fly-pg-may26-d2.100samples.decomposed.svs.rds`

### Post-decomposition/calling analysis

1. Cluster SVs into SV sites (for both SVs in the pangenome and calls across 100 samples): [`cluster-svs.R`](cluster-svs.R)
    - This will just add new columns to the calls to specify the assigned SV site.
    - Input: 
        - `<SAMPLE>.decomposed.svs.rds` for each sample in the pangenome
        - `fly-pg-may26-d2.100samples.decomposed.svs.rds`
    - Output:
        - `fly-pg-may26.svs.site.rol90.insd100.rds`
        - `fly-pg-may26-d2.100samples.decomposed.svs.site.rol90.insd100.rds`
    - The output files were deposited in:
        - [`http://public.gi.ucsc.edu/~jmonlong/hprc/minigraph-cactus-manuscript/fly/fly-pg-may26.svs.site.rol90.insd100.rds`](http://public.gi.ucsc.edu/~jmonlong/hprc/minigraph-cactus-manuscript/fly/fly-pg-may26.svs.site.rol90.insd100.rds)
        - [`http://public.gi.ucsc.edu/~jmonlong/hprc/minigraph-cactus-manuscript/fly/fly-pg-may26-d2.100samples.decomposed.svs.site.rol90.insd100.rds`](http://public.gi.ucsc.edu/~jmonlong/hprc/minigraph-cactus-manuscript/fly/fly-pg-may26-d2.100samples.decomposed.svs.site.rol90.insd100.rds)
2. Annotate the alleles with RepeatMasker: [`annotate-repeats.R`](annotate-repeats.R)
    - Annotates the major allele of each SV site
    - Uses the following docker container with RepeatMasker: `jmonlong/repeatmasker:release-4.0.9-p2`
    - Input: `fly-pg-may26-d2.100samples.decomposed.svs.site.rol90.insd100.rds`
    - Output: `fly-pg-may26-d2.100samples.decomposed.svs.site.rol90.insd100.rmsk.rds`
    - The output files were deposited at [`http://public.gi.ucsc.edu/~jmonlong/hprc/minigraph-cactus-manuscript/fly/fly-pg-may26-d2.100samples.decomposed.svs.site.rol90.insd100.rmsk.rds`](http://public.gi.ucsc.edu/~jmonlong/hprc/minigraph-cactus-manuscript/fly/fly-pg-may26-d2.100samples.decomposed.svs.site.rol90.insd100.rmsk.rds)
1. Make figures: [`svs-analysis.R`](svs-analysis.R).
    - Input:
        - `fly-pg-may26.svs.site.rol90.insd100.rds`
        - `fly-pg-may26-d2.100samples.decomposed.svs.site.rol90.insd100.rds`
        - `fly-pg-may26-d2.100samples.decomposed.svs.site.rol90.insd100.rmsk.rds`

