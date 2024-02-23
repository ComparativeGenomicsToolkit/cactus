# READ THIS FIRST

https://github.com/marbl/Primates#data-reuse-and-license

# 10-way T2T Ape Pangenome (human, chimp, bonobo, gorilla)

## Cactus Alignment

Cactus version: v2.7.1

Note that this uses the new `--noSplit` option that puts all contigs into one big graph without doing any kind of chromosome decomposition (as was done for the HPRC).  Also, we're requesting tons of memory (see log for actual usage).  Most processes don't need so much but the distance indexing takes over a terrabyte for this graph. 

In terms of graph construction, the 4 orang genomes can be added without much issue (I've tried this on the older assemblies). Not sure about the indexing though.  

```
TOIL_SLURM_ARGS="--partition=long --time=8000" cactus-pangenome ./js-pg ./10-t2t-apes-mc-2023v2.seqfile --outDir 10-t2t-apes-mc-2023v2 --outName 10-t2t-apes-mc-2023v2 --reference hs1 hg38 --noSplit --gbz clip full --gfa clip full --xg clip full --odgi --vcf --giraffe clip --haplo clip --vcfReference hs1 hg38 --logFile 10-t2t-apes-mc-2023v2.log  --batchSystem slurm --coordinationDir /data/tmp --caching false --batchLogsDir ./batch-logs --consMemory 1500Gi --indexMemory 1500Gi --mgMemory 500Gi --mgCores 72 --mapCores 8 --consCores 128 --indexCores 72 --giraffe clip
```

## VCF PostPorcessing

The same processing was used to make the "decomposed" HPRC vcfs.  Note that this process will filter out giant SVs as a 100kb threshold was used.

```
docker run -it --rm -v $(pwd)/:/data --user $(id -u):$(id -g) ghcr.io/pangenome/pggb:202402032147026ffe7f bash -c /data/vcf-bubwave.sh /data/10-t2t-apes-mc-2023v2.vcf.gz /data/10-t2t-apes-mc-2023v2-vcfbub.a100k.wave.vcf.gz
```

Left-align indels

```
vg paths -x 10-t2t-apes-mc-2023v2.xg -FQ hs1 | sed -e "s/hs1#0#//g" | bgzip > hs1.fa.gz
samtools faidx hs1.fa.gz
bcftools norm -f hs1.fa.gz 10-t2t-apes-mc-2023v2-vcfbub.a100k.wave.norm.vcf.gz | bcftools sort | bgzip > 10-t2t-apes-mc-2023v2-vcfbub.a100k.wave.norm.vcf.gz
tabix -fp vcf 10-t2t-apes-mc-2023v2-vcfbub.a100k.wave.norm.vcf.gz
```

## MAF Export for human

halRenameGenomes 10-t2t-apes-mc-2023v2.full.hal rename-hal-to-gca.tsv

for i in hs1 hg38 ; do TOIL_SLURM_ARGS="--partition=long --time=8000" cactus-hal2maf ./js ./10-t2t-apes-mc-2023v2.full.hal ./10-t2t-apes-mc-2023v2.${i}.maf.gz  --filterGapCausingDupes --refGenome $i --chunkSize 500000 --batchCores 64 --noAncestors --batchCount 16  --batchSystem slurm --caching false --logFile ./10-t2t-apes-mc-2023v2.${i}.maf.gz.log --batchLogsDir batch-logs-16apes --coordinationDir /data/tmp ;done

Then make a bigmaf (Note: using cactus commit a8bd77e65d7f7c26fd7a6d69a110d1fe23b275c9 for this one -- reproduce with v2.7.2)

for i in hs1 hg38; do TOIL_SLURM_ARGS="--partition=long --time=8000" cactus-maf2bigmaf ./js-bb ./10-t2t-apes-mc-2023v2.${i}.maf.gz ./10-t2t-apes-mc-2023v2.${i}.bigmaf.bb --refGenome $i --halFile ./10-t2t-apes-mc-2023v2.full.hal --logFile ./10-t2t-apes-mc-2023v2.${i}.bigmaf.bb.log  --batchLogsDir batch-logs-8apes --coordinationDir /data/tmp --batchSystem slurm; done


## Naming

The GCA accessions don't play nicely with the pangenome haplotype naming convetions, the alignment was generated with human-readable names.  To flip back to accessions in the HAL

```
halRenameGenomes 10-t2t-apes-mc-2023v2.full.hal rename-hal-to-gca.tsv
```

And to go back
```
halRenameGenomes 10-t2t-apes-mc-2023v2.full.hal rename-hal-from-gca.tsv
```

For the maf's, use

```
zcat 10-t2t-apes-mc-2023v2.hs1.maf.gz | sed -e's/hg002.1/GCA_018852605.2/g' \
-e 's/hg002.2/GCA_018852615.2/g' \
-e 's/mPanTro3.1/GCA_028858775.2/g' \
-e 's/mPanTro3.2/GCA_028858805.2/g' \
-e 's/mPanPan1.1/GCA_028858825.2/g' \
-e 's/mPanPan1.2/GCA_028858845.2/g' \
-e 's/mGorGor1.1/GCA_028885475.2/g' \
-e 's/mGorGor1.2/GCA_028885495.2/g' | bgzip > 10-t2t-apes-mc-2023v2.hs1.remamed.maf.gz
```

and

```
zcat 10-t2t-apes-mc-2023v2.hs1.renamed.maf.gz | sed -e's/GCA_018852605.2/hg002.1/g' \
-e 's/GCA_018852615.2/hg002.2/g' \
-e 's/GCA_028858775.2/mPanTro3.1/g' \
-e 's/GCA_028858805.2/mPanTro3.2/g' \
-e 's/GCA_028858825.2/mPanPan1.1/g' \
-e 's/GCA_028858845.2/mPanPan1.2/g' \
-e 's/GCA_028885475.2/mGorGor1.1/g' \
-e 's/GCA_028885495.2/mGorGor1.2/g' | bgzip > 10-t2t-apes-mc-2023v2.hs1.remamed.maf.gz
```
