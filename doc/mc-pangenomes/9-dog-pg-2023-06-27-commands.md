This is the same input data as "9-dog-pg-2022-09-23", just run with the current version of Cactus on SLURM.

Cactus version: v2.6.1

Note: `--consMemory 256Gi` probably very conservative and could be lowered (see memory usage in `9-dog-pg-2023-06-27.log`)

```
cactus-pangenome ./js ./9-dog-pg-2023-06-27-seqfile.txt --outName 9-dog-pg-2023-06-27 --outDir 9-dog-pg-2023-06-27 --reference canFam4 --batchSystem slurm --indexCores 64 --consCores 64 --mgCores 64 --mapCores 8 --logFile 9-dog-pg-2023-06-27.log --gbz --gfa --vcf --giraffe --chrom-og --chrom-vg --viz --consMemory 256Gi
```

Oops! (gaf file extension not added due to bug in this release)
```
mv 9-dog-pg-2023-06-27/9-dog-pg-2023-06-27 9-dog-pg-2023-06-27/9-dog-pg-2023-06-27.gaf.gz
```

