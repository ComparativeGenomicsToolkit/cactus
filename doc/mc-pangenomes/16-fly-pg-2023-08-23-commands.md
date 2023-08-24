# Steps to reproduce 16-fly-pg-2023-08-23

First, pull in the input assemblies that we used for the pangeome in the mc paper
```
wget https://s3-us-west-2.amazonaws.com/human-pangenomics/publications/mc_2022/mc_pangenomes/16-fruitfly-mc-2022-05-26/16-fly-softmasked-fa.tar.gz
tar zxf 16-fly-softmasked-fa.tar.gz
```

Make and index the pangenome with Cactus v2.6.7 on single machine (I ran on mustard.prism)
(Todo: this command does not work on slurm, because the 9Gb memory it estimates for minigraph construction is apparently too low)

Notable options:
* `--permissiveContigFilter` to help assign tiny contigs to reference chromosome
* `--haplo` to make the new subsampling index (.hapl) for giraffe, which obviates need for the filter (.d2) graph


```
cactus-pangenome ./js ./16-fly-pg-2023-08-23-seqfile.txt --outDir 16-fly-pg-2023-08-23 --outName 16-fly-pg-2023-08-23-seqfile.txt --reference dm6 --giraffe clip filter --gbz clip filter full --gfa clip filter full --vcf --permissiveContigFilter --haplo --chrom-vg clip filter --chrom-og full --viz --consCores 32 --indexCores 32 --mgCores 64 --mapCores 8 --logFile 16-fly-pg-2023-08-23.log  2> 16-fly-pg-2023-08-23.stderr
```

Oops, `--haplo` doesn't work with this data, apparently.  Until vg is patched, we run without

```
cactus-graphmap-join ./js --vg 16-fly-pg-2023-08-23/chrom-alignments/*.vg --hal 16-fly-pg-2023-08-23/chrom-alignments/*.hal --outDir 16-fly-pg-2023-08-23 --outName 16-fly-pg-2023-08-23-seqfile.txt --reference dm6 --giraffe clip filter --gbz clip filter full --gfa clip filter full --vcf --chrom-vg clip filter --chrom-og full --viz --indexCores 32 --logFile 16-fly-pg-2023-08-23.join.log  2> 16-fly-pg-2023-08-23.join.stderr
```
