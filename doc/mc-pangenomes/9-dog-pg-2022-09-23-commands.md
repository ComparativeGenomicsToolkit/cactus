```

export MYBUCKET=s3://vg-k8s/users/hickey/dog-pangenome/
export MYJOBSTORE=aws:us-west-2:cactus-hprc-jobstore-dog
export VERSION=sep23

cactus-minigraph ${MYJOBSTORE} 9-dog-pg.txt ${MYBUCKET}/dog-pg-${VERSION}.gfa.gz --reference canFam4 --logFile dog-pg-${VERSION}-minigraph.log --batchSystem mesos --provisioner aws --defaultPreemptable --nodeType r5.8xlarge:1.5 --nodeStorage 650 --maxNodes 25 --betaInertia 0 --targetTime 1 --mapCores 32

cactus-graphmap ${MYJOBSTORE} 9-dog-pg.txt ${MYBUCKET}/dog-pg-${VERSION}.gfa.gz ${MYBUCKET}/dog-pg-${VERSION}.paf  --outputFasta ${MYBUCKET}/dog-pg-${VERSION}.gfa.fa.gz --reference canFam4 --mapCores 16 --delFilter 10000000 --logFile dog-pg-${VERSION}-graphmap.log  --batchSystem mesos --provisioner aws --defaultPreemptable --nodeType r5.8xlarge:1.5 --nodeStorage 650 --maxNodes 25 --betaInertia 0 --targetTime 1 --mapCores 16

cactus-graphmap-split ${MYJOBSTORE} 9-dog-pg.txt  ${MYBUCKET}/dog-pg-${VERSION}.gfa.gz ${MYBUCKET}/dog-pg-${VERSION}.paf --outDir ${MYBUCKET}/chroms-dog-pg-${VERSION} --otherContig chrOther --refContigs $(for i in `seq 38`; do echo chr$i; done ; echo "chrX chrM") --reference canFam4 --logFile dog-pg-${VERSION}-split.log  --batchSystem mesos --provisioner aws --defaultPreemptable --nodeType r5.8xlarge:1.5 --nodeStorage 650 --maxNodes 25 --betaInertia 0 --targetTime 1

aws s3 cp ${MYBUCKET}/chroms-dog-pg-${VERSION}/chromfile.txt .
cactus-align-batch ${MYJOBSTORE} ./chromfile.txt ${MYBUCKET}/align-dog-pg-${VERSION} --alignCores 16 --realTimeLogging --alignOptions "--pangenome --maxLen 10000 --reference canFam4 --outVG" --logFile dog-pg-${VERSION}-align.log --batchSystem mesos --provisioner aws --defaultPreemptable --nodeType r5.8xlarge:1.5 --nodeStorage 1000 --maxNodes 20 --betaInertia 0 --targetTime 1

cactus-graphmap-join ${MYJOBSTORE} --vg $(for j in $(for i in `seq 38`; do echo chr$i; done ; echo "chrX chrM chrOther"); do echo ${MYBUCKET}/align-dog-pg-${VERSION}/${j}.vg; done) --hal $(for j in $(for i in `seq 38`; do echo chr$i; done ; echo "chrX chrM chrOther"); do echo ${MYBUCKET}/align-dog-pg-${VERSION}/${j}.hal; done) --outDir ${MYBUCKET}/ --outName dog-pg-${VERSION}-full --reference canFam4 --gfaffix  --wlineSep "."  --logFile dog-pg-${VERSION}-join-full.log --batchSystem mesos --provisioner aws --defaultPreemptable --nodeType r5.8xlarge:1.5 --nodeStorage 1000 --maxNodes 1 --indexCores 31 --realTimeLogging

cactus-graphmap-join ${MYJOBSTORE} --vg $(for j in $(for i in `seq 38`; do echo chr$i; done ; echo "chrX chrM chrOther"); do echo ${MYBUCKET}/clip-dog-pg-${VERSION}-full/${j}.vg; done) --outDir ${MYBUCKET}/ --outName dog-pg-${VERSION} --reference canFam4  --wlineSep "." --clipLength 10000 --clipNonMinigraph --vcf --giraffe --preserveIDs --logFile dog-pg-${VERSION}-join.log --batchSystem mesos --provisioner aws --defaultPreemptable --nodeType r5.8xlarge:1.5 --nodeStorage 1000 --maxNodes 2 --indexCores 32 --betaInertia 0 --targetTime 1 --realTimeLogging

cactus-graphmap-join ${MYJOBSTORE} --vg $(for j in $(for i in `seq 38`; do echo chr$i; done ; echo "chrX chrM chrOther"); do echo ${MYBUCKET}/clip-dog-pg-${VERSION}/${j}.vg; done) --outDir ${MYBUCKET} --outName dog-pg-${VERSION}-minaf.0.22 --reference canFam4  --wlineSep "." --vgClipOpts "-d 2 -m 1000" --preserveIDs --giraffe --logFile dog-pg-${VERSION}-join-minaf.0.22.log --batchSystem mesos --provisioner aws --defaultPreemptable --nodeType r5.8xlarge:1.5 --nodeStorage 1000 --maxNodes 2 --indexCores 32 --realTimeLogging

for i in *.log ; do aws s3 cp $i ${MYBUCKET}/doc/$i; done

```
