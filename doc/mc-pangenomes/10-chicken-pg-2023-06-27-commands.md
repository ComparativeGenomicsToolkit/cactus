This is the same input data as "10-chicken-pg-2022-09-23", just run with the current version of Cactus on SLURM.

Cactus version: v2.6.1

Note: `--consMemory 256Gi` is very conservative and could be lowered by at least half(see memory usage in `10-chicken-pg-2023-06-27.log`).  But unfortunately `--consMemory` is required to run this data on slurm, as some of the smaller chromosomes, ex chr30, use more memory than Cactus estimates based on their tiny size.  If not running on slurm, this won't be an issue.

```
cactus-pangenome ./js ./10-chicken-pg-2022-09-23-seqfile.txt --outName 10-chicken-pg-2023-06-7 --outDir 10-chicken-pg-2023-06-27 --reference galGal6 --batchSystem slurm --indexCores 64 --consCores 64 --mgCores 64 --mapCores 8 --logFile 10-chicken-pg-2023-06-27.log --maxLocalJobs 1000 --gbz --gfa --vcf --giraffe --chrom-og --chrom-vg --viz --consMemory 256Gi
```

Oops! (export of debug file broken in this release)
```
mv 10-chicken-pg-2023-06-27/10-chicken-pg-2023-06-27 10-chicken-pg-2023-06-27/10-chicken-pg-2023-06-27.paf.unfiltered.gz
```
