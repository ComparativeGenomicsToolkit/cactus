```

export MYBUCKET=s3://vg-k8s/users/hickey/soybean-pangenome
export MYJOBSTORE=aws:us-west-2:cactus-hprc-jobstore-soybean
export VERSION=sep26

wget https://hgdownload.soe.ucsc.edu/hubs/GCF/000/004/515/GCF_000004515.6/GCF_000004515.6.fa.gz
wget https://hgdownload.soe.ucsc.edu/hubs/GCF/000/004/515/GCF_000004515.6/GCF_000004515.6.chromAlias.txt
gzip -d GCF_000004515.6.fa.gz
awk 'FNR==NR{a[">"$1]=$(NF);next}{$1=($1 in a)?">"a[$1]:$1}1' GCF_000004515.6.chromAlias.txt GCF_000004515.6.fa | bgzip --threads 8 > GCF_000004515.6.ucsc.fa.gz
aws s3 cp GCF_000004515.6.ucsc.fa.gz $MYBUCKET/
sed -i  17-soybean-pg.txt  -e "s#https://hgdownload.soe.ucsc.edu/hubs/GCF/000/004/515/GCF_000004515.6/GCF_000004515.6.fa.gz#${MYBUCKET}/GCF_000004515.6.ucsc.fa.gz#g"

cactus-minigraph ${MYJOBSTORE} 17-soybean-pg.txt ${MYBUCKET}/soybean-pg-${VERSION}.gfa.gz --reference Glycine_max_v4_0 --logFile soybean-pg-${VERSION}-minigraph.log --batchSystem mesos --provisioner aws --defaultPreemptable --nodeType r5.8xlarge:1.5 --nodeStorage 650 --maxNodes 25 --betaInertia 0 --targetTime 1 --mapCores 32

cactus-graphmap ${MYJOBSTORE} 17-soybean-pg.txt ${MYBUCKET}/soybean-pg-${VERSION}.gfa.gz ${MYBUCKET}/soybean-pg-${VERSION}.paf  --outputFasta ${MYBUCKET}/soybean-pg-${VERSION}.gfa.fa.gz --reference Glycine_max_v4_0 --mapCores 16 --delFilter 10000000 --logFile soybean-pg-${VERSION}-graphmap.log  --batchSystem mesos --provisioner aws --defaultPreemptable --nodeType r5.8xlarge:1.5 --nodeStorage 650 --maxNodes 25 --betaInertia 0 --targetTime 1 --mapCores 16

cactus-graphmap-split ${MYJOBSTORE} 17-soybean-pg.txt  ${MYBUCKET}/soybean-pg-${VERSION}.gfa.gz ${MYBUCKET}/soybean-pg-${VERSION}.paf --outDir ${MYBUCKET}/chroms-soybean-pg-${VERSION} --otherContig chrOther --refContigs $(for i in `seq 20`; do echo chr$i; done ; echo "chrM") --reference Glycine_max_v4_0 --logFile soybean-pg-${VERSION}-split.log  --batchSystem mesos --provisioner aws --defaultPreemptable --nodeType r5.8xlarge:1.5 --nodeStorage 650 --maxNodes 25 --betaInertia 0 --targetTime 1

aws s3 cp ${MYBUCKET}/chroms-soybean-pg-${VERSION}/chromfile.txt .
cactus-align-batch ${MYJOBSTORE} ./chromfile.txt ${MYBUCKET}/align-soybean-pg-${VERSION} --alignCores 16 --realTimeLogging --alignOptions "--pangenome --maxLen 10000 --reference Glycine_max_v4_0 --outVG" --logFile soybean-pg-${VERSION}-align.log --batchSystem mesos --provisioner aws --defaultPreemptable --nodeType r5.8xlarge:1.5 --nodeStorage 1000 --maxNodes 20 --betaInertia 0 --targetTime 1

cactus-graphmap-join ${MYJOBSTORE} --vg $(for j in $(for i in `seq 20`; do echo chr$i; done ; echo "chrM chrOther"); do echo ${MYBUCKET}/align-soybean-pg-${VERSION}/${j}.vg; done) --hal $(for j in $(for i in `seq 20`; do echo chr$i; done ; echo "chrM chrOther"); do echo ${MYBUCKET}/align-soybean-pg-${VERSION}/${j}.hal; done) --outDir ${MYBUCKET}/ --outName soybean-pg-${VERSION}-full --reference Glycine_max_v4_0 --gfaffix  --wlineSep "."  --logFile soybean-pg-${VERSION}-join-full.log --batchSystem mesos --provisioner aws --defaultPreemptable --nodeType r5.16xlarge --nodeStorage 1000 --maxNodes 1 --indexCores 63 --realTimeLogging

cactus-graphmap-join ${MYJOBSTORE} --vg $(for j in $(for i in `seq 20`; do echo chr$i; done ; echo "chrM chrOther"); do echo ${MYBUCKET}/clip-soybean-pg-${VERSION}-full/${j}.vg; done) --outDir ${MYBUCKET}/ --outName soybean-pg-${VERSION} --reference Glycine_max_v4_0  --wlineSep "." --clipLength 10000 --clipNonMinigraph --vcf --giraffe --preserveIDs --logFile soybean-pg-${VERSION}-join.log --batchSystem mesos --provisioner aws --defaultPreemptable --nodeType r5.8xlarge --nodeStorage 1000 --maxNodes 2 --indexCores 32 --betaInertia 0 --targetTime 1 --realTimeLogging

cactus-graphmap-join ${MYJOBSTORE} --vg $(for j in $(for i in `seq 20`; do echo chr$i; done ; echo "chrM chrOther"); do echo ${MYBUCKET}/clip-soybean-pg-${VERSION}/${j}.vg; done) --outDir ${MYBUCKET} --outName soybean-pg-${VERSION}-minaf.0.22 --reference Glycine_max_v4_0  --wlineSep "." --vgClipOpts "-d 2 -m 1000" --preserveIDs --giraffe --logFile soybean-pg-${VERSION}-minaf.0.22.log --batchSystem mesos --provisioner aws --defaultPreemptable --nodeType r5.16xlarge --nodeStorage 1000 --maxNodes 2 --indexCores 64 --realTimeLogging

for i in *.log ; do aws s3 cp $i ${MYBUCKET}/doc/$i; done

```
