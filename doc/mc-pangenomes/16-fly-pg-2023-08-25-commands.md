# Steps to reproduce 16-fly-pg-2023-08-25

First, pull in the input assemblies that we used for the pangeome in the mc paper
```
wget https://s3-us-west-2.amazonaws.com/human-pangenomics/publications/mc_2022/mc_pangenomes/16-fruitfly-mc-2022-05-26/16-fly-softmasked-fa.tar.gz
tar zxf 16-fly-softmasked-fa.tar.gz
```

Make and index the pangenome with Cactus commit 1f537de69a9b7cb2c6f46ecc048023f92f911ed2 on slurm (which contains some fixes that will appear in v2.6.8)  

Notable options:
* `--permissiveContigFilter` to help assign tiny contigs to reference chromosome
* `--haplo` to make the new subsampling index (.hapl) for giraffe, which obviates need for the filter (.d2) graph


```
cactus-pangenome ./js ./16-fly-pg-2023-08-25-seqfile.txt --outDir 16-fly-pg-2023-08-25 --outName 16-fly-pg-2023-08-25 --reference dm6 --giraffe clip filter --gbz clip filter full --gfa clip filter full --vcf --permissiveContigFilter --haplo --chrom-vg clip filter --chrom-og full --viz --consCores 32 --indexCores 32 --mgCores 64 --mapCores 8 --batchSystem slurm --logFile 16-fly-pg-2023-08-25.log  2> 16-fly-pg-2023-08-25.stderr
```

