# 4-way T2T Orangutan Pangenome (mPonAbe1, mPonPyg2)

## Cactus Alignment

Cactus version: v2.8.1

Note that this uses the new `--noSplit` option that puts all contigs into one big graph without doing any kind of chromosome decomposition (as was done for the HPRC).  Also, we're requesting tons of memory (see log for actual usage).

TODO: Can we run chromosomes independently for inter-orang alignments?


```
TOIL_SLURM_ARGS="--partition=long --time=8000" cactus-pangenome ./js-pg ./4-t2t-orangs-mc-2023v2.seqfile --outDir 4-t2t-orangs-mc-2023v2 --outName 4-t2t-orangs-mc-2023v2 --reference mPonAbe1_pri mPonAbe1_alt --noSplit --gbz clip full --gfa clip full --xg clip full --odgi --vcf --giraffe clip --haplo clip --vcfReference mPonAbe1_pri mPonAbe1_alt --logFile 4-t2t-orangs-mc-2023v2.log  --batchSystem slurm --coordinationDir /data/tmp --batchLogsDir ./batch-logs --consMemory 1500Gi --indexMemory 1500Gi --mgMemory 500Gi --mgCores 72 --mapCores 8 --consCores 128 --indexCores 72 --giraffe clip
```

## VCF PostPorcessing

The same processing was used to make the "decomposed" HPRC vcfs.  Note that this process will filter out giant SVs as a 100kb threshold was used.

```
docker run -it --rm -v $(pwd)/:/data --user $(id -u):$(id -g) ghcr.io/pangenome/pggb:202402032147026ffe7f bash -c "/data/vcf-bubwave.sh /data/4-t2t-orangs-mc-2023v2.raw.vcf.gz /data/4-t2t-orangs-mc-2023v2-vcfbub.a100k.wave.vcf.gz"
```

Left-align indels

```
vg paths -x 4-t2t-orangs-mc-2023v2.xg -FQ mPonAbe1_pri | sed -e "s/mPonAbe1_pri#0#//g" | bgzip > mPonAbe1_pri.fa.gz
samtools faidx mPonAbe1_pri.fa.gz
bcftools norm -f mPonAbe1_pri.fa.gz 4-t2t-orangs-mc-2023v2-vcfbub.a100k.wave.vcf.gz | bcftools sort | bgzip > 4-t2t-orangs-mc-2023v2-vcfbub.a100k.wave.norm.vcf.gz
tabix -fp vcf 4-t2t-orangs-mc-2023v2-vcfbub.a100k.wave.norm.vcf.gz
```

## MAF

```
halRenameGenomes 4-t2t-orangs-mc-2023v2.full.hal rename-hal-to-gca.tsv
```
```
for i in GCA_028885655.2 GCA_028885685.2 ; do TOIL_SLURM_ARGS="--partition=long --time=8000" cactus-hal2maf ./js ./4-t2t-orangs-mc-2023v2.full.hal ./4-t2t-orangs-mc-2023v2.${i}.maf.gz  --filterGapCausingDupes --refGenome $i --chunkSize 500000 --batchCores 64 --noAncestors --batchCount 16  --batchSystem slurm --caching false --logFile ./4-t2t-orangs-mc-2023v2.${i}.maf.gz.log --batchLogsDir batch-logs-16apes --coordinationDir /data/tmp ; done
```

Then make a bigmaf 
```
for i in hs1 hg38; do TOIL_SLURM_ARGS="--partition=long --time=8000" cactus-maf2bigmaf ./js-bb ./4-t2t-orangs-mc-2023v2.${i}.maf.gz ./10-t2t-apes-mc-2023v2.${i}.bigmaf.bb --refGenome $i --halFile ./4-t2t-orangs-mc-2023v2.full.hal --logFile ./4-t2t-orangs-mc-2023v2.${i}.bigmaf.bb.log  --batchLogsDir batch-logs-8apes --coordinationDir /data/tmp --batchSystem slurm; done
```

## Chains

```
TOIL_SLURM_ARGS="--partition=long --time=8000" cactus-hal2chains ./js_chains10 ./4-t2t-orangs-mc-2023v2.full.hal ./4-t2t-orangs-mc-2023v2-chains  --inMemory --bigChain  --caching false --logFile ./4-t2t-orangs-mc-2023v2-chains.log --batchLogsDir batch-logs --batchSystem slurm --coordinationDir /data/tmp
```

```
halRenameGenomes 4-t2t-orangs-mc-2023v2.full.hal rename-hal-from-gca.tsv
```
