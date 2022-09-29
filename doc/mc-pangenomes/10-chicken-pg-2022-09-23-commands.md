```

export MYBUCKET=s3://vg-k8s/users/hickey/chicken-pangenome/
export MYJOBSTORE=aws:us-west-2:cactus-hprc-jobstore-chicken
export VERSION=sep23

cactus-minigraph ${MYJOBSTORE} 10-chicken-pg.txt ${MYBUCKET}/chicken-pg-${VERSION}.gfa.gz --reference galGal6 --logFile chicken-pg-${VERSION}-minigraph.log --batchSystem mesos --provisioner aws --defaultPreemptable --nodeType r5.8xlarge:1.5 --nodeStorage 650 --maxNodes 25 --betaInertia 0 --targetTime 1 --mapCores 32

cactus-graphmap ${MYJOBSTORE} 10-chicken-pg.txt ${MYBUCKET}/chicken-pg-${VERSION}.gfa.gz ${MYBUCKET}/chicken-pg-${VERSION}.paf  --outputFasta ${MYBUCKET}/chicken-pg-${VERSION}.gfa.fa.gz --reference galGal6 --mapCores 16 --delFilter 10000000 --logFile chicken-pg-${VERSION}-graphmap.log  --batchSystem mesos --provisioner aws --defaultPreemptable --nodeType r5.8xlarge:1.5 --nodeStorage 650 --maxNodes 25 --betaInertia 0 --targetTime 1 --mapCores 16

cactus-graphmap-split ${MYJOBSTORE} 10-chicken-pg.txt  ${MYBUCKET}/chicken-pg-${VERSION}.gfa.gz ${MYBUCKET}/chicken-pg-${VERSION}.paf --outDir ${MYBUCKET}/chroms-chicken-pg-${VERSION} --otherContig chrOther --refContigs $(for i in `seq 33`; do echo chr$i; done ; echo "chrW chrM") --reference galGal6 --logFile chicken-pg-${VERSION}-split.log  --batchSystem mesos --provisioner aws --defaultPreemptable --nodeType r5.8xlarge:1.5 --nodeStorage 650 --maxNodes 25 --betaInertia 0 --targetTime 1

aws s3 cp ${MYBUCKET}/chroms-chicken-pg-${VERSION}/chromfile.txt .
cactus-align-batch ${MYJOBSTORE} ./chromfile.txt ${MYBUCKET}/align-chicken-pg-${VERSION} --alignCores 16 --realTimeLogging --alignOptions "--pangenome --maxLen 10000 --reference galGal6 --outVG" --logFile chicken-pg-${VERSION}-align.log --batchSystem mesos --provisioner aws --defaultPreemptable --nodeType r5.8xlarge:1.5 --nodeStorage 1000 --maxNodes 20 --betaInertia 0 --targetTime 1

cactus-graphmap-join ${MYJOBSTORE} --vg $(for j in $(for i in `seq 28`; do echo chr$i; done ; echo "chr30 chr31 chr32 chr33 chrW chrM chrOther"); do echo ${MYBUCKET}/align-chicken-pg-${VERSION}/${j}.vg; done) --hal $(for j in $(for i in `seq 28`; do echo chr$i; done ; echo "chr30 chr31 chr32 chr33 chrW chrM chrOther"); do echo ${MYBUCKET}/align-chicken-pg-${VERSION}/${j}.hal; done) --outDir ${MYBUCKET}/ --outName chicken-pg-${VERSION}-full --reference galGal6 --gfaffix  --wlineSep "."  --logFile chicken-pg-${VERSION}-join-full.log --batchSystem mesos --provisioner aws --defaultPreemptable --nodeType r5.8xlarge:1.5 --nodeStorage 1000 --maxNodes 1 --indexCores 31 --realTimeLogging

cactus-graphmap-join ${MYJOBSTORE} --vg $(for j in $(for i in `seq 28`; do echo chr$i; done ; echo "chr30 chr31 chr32 chr33 chrW chrM chrOther"); do echo ${MYBUCKET}/clip-chicken-pg-${VERSION}-full/${j}.vg; done) --outDir ${MYBUCKET}/ --outName chicken-pg-${VERSION} --reference galGal6  --wlineSep "." --clipLength 10000 --clipNonMinigraph --vcf --giraffe --preserveIDs --logFile chicken-pg-${VERSION}-join.log --batchSystem mesos --provisioner aws --defaultPreemptable --nodeType r5.8xlarge:1.5 --nodeStorage 1000 --maxNodes 2 --indexCores 32 --betaInertia 0 --targetTime 1 --realTimeLogging

cactus-graphmap-join ${MYJOBSTORE} --vg $(for j in $(for i in `seq 28`; do echo chr$i; done ; echo "chr30 chr31 chr32 chr33 chrW chrM chrOther"); do echo ${MYBUCKET}/clip-chicken-pg-${VERSION}/${j}.vg; done) --outDir ${MYBUCKET} --outName chicken-pg-${VERSION}-minaf.0.2 --reference galGal6  --wlineSep "." --vgClipOpts "-d 2 -m 1000" --preserveIDs --giraffe --logFile chicken-pg-${VERSION}-minaf.0.2.log --batchSystem mesos --provisioner aws --defaultPreemptable --nodeType r5.8xlarge:1.5 --nodeStorage 1000 --maxNodes 2 --indexCores 32 --realTimeLogging 

for i in *.log ; do aws s3 cp $i ${MYBUCKET}/doc/$i; done

```
