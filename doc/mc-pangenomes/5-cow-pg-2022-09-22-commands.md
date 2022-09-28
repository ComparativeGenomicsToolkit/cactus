```

export MYBUCKET=s3://vg-k8s/users/hickey/cow-pangenome/
export MYJOBSTORE=aws:us-west-2:cactus-hprc-jobstore-cow
export VERSION=sep22

cactus-minigraph ${MYJOBSTORE} 5-cow-pg.txt ${MYBUCKET}/cow-pg-${VERSION}.gfa.gz --reference bosTau9 --logFile cow-pg-${VERSION}-minigraph.log --batchSystem mesos --provisioner aws --defaultPreemptable --nodeType r5.8xlarge:1.5 --nodeStorage 650 --maxNodes 25 --betaInertia 0 --targetTime 1 --mapCores 16

cactus-graphmap ${MYJOBSTORE} 5-cow-pg.txt ${MYBUCKET}/cow-pg-${VERSION}.gfa.gz ${MYBUCKET}/cow-pg-${VERSION}.paf  --outputFasta ${MYBUCKET}/cow-pg-${VERSION}.gfa.fa.gz --reference bosTau9 --mapCores 16 --delFilter 10000000 --logFile cow-pg-${VERSION}-graphmap.log  --batchSystem mesos --provisioner aws --defaultPreemptable --nodeType r5.8xlarge:1.5 --nodeStorage 650 --maxNodes 25 --betaInertia 0 --targetTime 1 --mapCores 16

cactus-graphmap-split ${MYJOBSTORE} 5-cow-pg.txt  ${MYBUCKET}/cow-pg-${VERSION}.gfa.gz ${MYBUCKET}/cow-pg-${VERSION}.paf --outDir ${MYBUCKET}/chroms-cow-pg-${VERSION} --otherContig chrOther --refContigs $(for i in `seq 29`; do echo chr$i; done ; echo "chrX chrM") --reference bosTau9 --logFile cow-pg-${VERSION}-split.log  --batchSystem mesos --provisioner aws --defaultPreemptable --nodeType r5.8xlarge:1.5 --nodeStorage 650 --maxNodes 25 --betaInertia 0 --targetTime 1

aws s3 cp ${MYBUCKET}/chroms-cow-pg-${VERSION}/chromfile.txt .
cactus-align-batch ${MYJOBSTORE} ./chromfile.txt ${MYBUCKET}/align-cow-pg-${VERSION} --alignCores 16 --realTimeLogging --alignOptions "--pangenome --maxLen 10000 --reference bosTau9 --outVG" --logFile cow-pg-${VERSION}-align.log --batchSystem mesos --provisioner aws --defaultPreemptable --nodeType r5.8xlarge:1.5 --nodeStorage 1000 --maxNodes 20 --betaInertia 0 --targetTime 1

cactus-graphmap-join ${MYJOBSTORE} --vg $(for j in $(for i in `seq 29`; do echo chr$i; done ; echo "chrX chrM chrOther"); do echo ${MYBUCKET}/align-cow-pg-${VERSION}/${j}.vg; done) --hal $(for j in $(for i in `seq 29`; do echo chr$i; done ; echo "chrX chrM chrOther"); do echo ${MYBUCKET}/align-cow-pg-${VERSION}/${j}.hal; done) --outDir ${MYBUCKET}/ --outName cow-pg-${VERSION}-full --reference bosTau9 --gfaffix  --wlineSep "."  --logFile cow-pg-${VERSION}-join-full.log --batchSystem mesos --provisioner aws --defaultPreemptable --nodeType r5.16xlarge --nodeStorage 1000 --maxNodes 1 --indexCores 63 --realTimeLogging

cactus-graphmap-join ${MYJOBSTORE} --vg $(for j in $(for i in `seq 29`; do echo chr$i; done ; echo "chrX chrM chrOther"); do echo ${MYBUCKET}/clip-cow-pg-${VERSION}-full/${j}.vg; done) --outDir ${MYBUCKET}/ --outName cow-pg-${VERSION} --reference bosTau9  --wlineSep "." --clipLength 10000 --clipNonMinigraph --vcf --giraffe --preserveIDs --logFile cow-pg-${VERSION}-join.log --batchSystem mesos --provisioner aws --defaultPreemptable --nodeType r5.8xlarge --nodeStorage 1000 --maxNodes 2 --indexCores 32 --betaInertia 0 --targetTime 1 --realTimeLogging

cactus-graphmap-join ${MYJOBSTORE} --vg $(for j in $(for i in `seq 29`; do echo chr$i; done ; echo "chrX chrM chrOther"); do echo ${MYBUCKET}/clip-cow-pg-${VERSION}/${j}.vg; done) --outDir ${MYBUCKET} --outName cow-pg-${VERSION}-minaf.0.4 --reference bosTau9  --wlineSep "." --vgClipOpts "-d 2 -m 1000" --preserveIDs --giraffe --logFile cow-pg-${VERSION}-minaf.0.4.log --batchSystem mesos --provisioner aws --defaultPreemptable --nodeType r5.16xlarge --nodeStorage 1000 --maxNodes 2 --indexCores 64 --realTimeLogging

for i in *.log ; do aws s3 cp $i ${MYBUCKET}/doc/$i; done


```
