```

export MYBUCKET=s3://vg-k8s/users/hickey/mouse-pangenome
export MYJOBSTORE=aws:us-west-2:cactus-hprc-jobstore-mouse
export VERSION=sep23

cactus-minigraph ${MYJOBSTORE} 30-mouse-pg.txt ${MYBUCKET}/mouse-pg-${VERSION}.gfa.gz --reference mm39 --logFile mouse-pg-${VERSION}-minigraph.log --batchSystem mesos --provisioner aws --defaultPreemptable --nodeType r5.4xlarge --nodeStorage 650 --maxNodes 25 --betaInertia 0 --targetTime 1 --mapCores 16

cactus-graphmap ${MYJOBSTORE} 30-mouse-pg.txt ${MYBUCKET}/mouse-pg-${VERSION}.gfa.gz ${MYBUCKET}/mouse-pg-${VERSION}.paf  --outputFasta ${MYBUCKET}/mouse-pg-${VERSION}.gfa.fa.gz --reference mm39 --mapCores 16 --delFilter 10000000 --logFile mouse-pg-${VERSION}-graphmap.log  --batchSystem mesos --provisioner aws --defaultPreemptable --nodeType r5.8xlarge:1.5 --nodeStorage 650 --maxNodes 25 --betaInertia 0 --targetTime 1 --mapCores 16

cactus-graphmap-split ${MYJOBSTORE} 30-mouse-pg.txt  ${MYBUCKET}/mouse-pg-${VERSION}.gfa.gz ${MYBUCKET}/mouse-pg-${VERSION}.paf --outDir ${MYBUCKET}/chroms-mouse-pg-${VERSION} --otherContig chrOther --refContigs $(for i in `seq 19`; do echo chr$i; done ; echo "chrX chrY chrM") --reference mm39 --logFile mouse-pg-${VERSION}-split.log  --batchSystem mesos --provisioner aws --defaultPreemptable --nodeType r5.8xlarge:1.5 --nodeStorage 650 --maxNodes 25 --betaInertia 0 --targetTime 1

aws s3 cp ${MYBUCKET}/chroms-mouse-pg-${VERSION}/chromfile.txt .
cactus-align-batch ${MYJOBSTORE} ./chromfile.txt ${MYBUCKET}/align-mouse-pg-${VERSION} --alignCores 16 --realTimeLogging --alignOptions "--pangenome --maxLen 10000 --reference mm39 --outVG" --logFile mouse-pg-${VERSION}-align.log --batchSystem mesos --provisioner aws --defaultPreemptable --nodeType r5.8xlarge:1.5 --nodeStorage 1000 --maxNodes 20 --betaInertia 0 --targetTime 1

cactus-graphmap-join ${MYJOBSTORE} --vg $(for j in $(for i in `seq 19`; do echo chr$i; done ; echo "chrX chrY chrM chrOther"); do echo ${MYBUCKET}/align-mouse-pg-${VERSION}/${j}.vg; done) --hal $(for j in $(for i in `seq 19`; do echo chr$i; done ; echo "chrX chrY chrM chrOther"); do echo ${MYBUCKET}/align-mouse-pg-${VERSION}/${j}.hal; done) --outDir ${MYBUCKET}/ --outName mouse-pg-${VERSION}-full --reference mm39 --gfaffix  --wlineSep "."  --logFile mouse-pg-${VERSION}-join-full.log --batchSystem mesos --provisioner aws --defaultPreemptable --nodeType r5.8xlarge:1.5 --nodeStorage 1000 --maxNodes 1 --indexCores 31 --realTimeLogging

cactus-graphmap-join ${MYJOBSTORE} --vg $(for j in $(for i in `seq 19`; do echo chr$i; done ; echo "chrX chrY chrM chrOther"); do echo ${MYBUCKET}/clip-mouse-pg-${VERSION}-full/${j}.vg; done) --outDir ${MYBUCKET}/ --outName mouse-pg-${VERSION} --reference mm39  --wlineSep "." --clipLength 10000 --clipNonMinigraph --vcf --giraffe --preserveIDs --logFile mouse-pg-${VERSION}-join.log --batchSystem mesos --provisioner aws --defaultPreemptable --nodeType r5.8xlarge:1.5 --nodeStorage 1000 --maxNodes 2 --indexCores 32 --betaInertia 0 --targetTime 1 --realTimeLogging

cactus-graphmap-join ${MYJOBSTORE} --vg $(for j in $(for i in `seq 19`; do echo chr$i; done ; echo "chrX chrY chrM chrOther"); do echo ${MYBUCKET}/clip-mouse-pg-${VERSION}/${j}.vg; done) --outDir ${MYBUCKET} --outName mouse-pg-${VERSION}-mc-grch19-minaf.0.10 --reference mm39  --wlineSep "." --vgClipOpts "-d 3 -m 1000" --preserveIDs --giraffe --logFile mouse-pg-${VERSION}-minaf.0.30.log --batchSystem mesos --provisioner aws --defaultPreemptable --nodeType r5.8xlarge:1.5 --nodeStorage 1000 --maxNodes 2 --indexCores 32 --realTimeLogging 

for i in *.log ; do aws s3 cp $i ${MYBUCKET}/doc/$i; done

```
