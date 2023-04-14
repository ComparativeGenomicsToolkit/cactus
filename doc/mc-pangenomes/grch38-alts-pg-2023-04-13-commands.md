Cactus version used: v2.5.0

Download hg38 and split it into one fasta file for non-alt contigs, and a seperate fast file for each alt.

```
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
gzip -d hg38.fa.gz
samtools faidx hg38.fa
mkdir -p fa
rm -f fa/hg38_no_alts.fa
for contig in $(grep -v _alt hg38.fa.fai | awk '{print $1}'); do samtools faidx hg38.fa $contig >> fa/hg38_no_alts.fa; done
bgzip fa/hg38_no_alts.fa --threads 16
printf "GRCh38\tfa/hg38_no_alts.fa.gz\n" > grch38_alts.seqfile
for contig in $(grep  _alt hg38.fa.fai | awk '{print $1}'); \
do samtools faidx hg38.fa $contig | bgzip > fa/${contig}.fa.gz; \
printf "${contig}\tfa/${contig}.fa.gz\n" >> grch38_alts.seqfile ; done
```

Make the pangenome, setting a very permissive overlap threshold to make sure all contigs get into a chromosome

```
mkdir -p work
cactus-pangenome ./js ./grch38_alts.seqfile --reference GRCh38 --gbz clip --giraffe clip --vcf --outName grch38-alts-apr13 --outDir grch38-alts --logFile grch38-alts-apr13.log --indexCores 32 --mapCores 8 --consCores 8 --refContigs $(for i in `seq 22`; do printf "chr$i "; done ; echo "chrX chrY chrM") --otherContig chrOther --permissiveContigFilter 0.05 --workDir work
```

Note: the actual command run had quotes around the refContigs, causing them all to be run at once. This slows it down but doesn't affect the output. Quotes are removed above so others don't make the same mistake.
