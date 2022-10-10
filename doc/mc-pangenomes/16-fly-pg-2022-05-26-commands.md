```
aws s3 cp s3://vg-k8s/users/hickey/fly-pangenome/16-fly-softmasked-fa.tar.gz .
tar zxf 16-fly-softmasked-fa.tar.gz

export VERSION=2022-05-26

cactus-minigraph ./js ./softmask.fly.txt 16-fly-pg-${VERSION}-minigraph.gfa.gz --reference dm6 --maxCores 32 --mapCores 32 --realTimeLogging --logFile 16-fly-pg-${VERSION}.minigraph.log

cactus-graphmap ./js ./softmask.fly.txt ./16-fly-pg-${VERSION}-minigraph.gfa.gz ./16-fly-pg-${VERSION}.paf --reference dm6 --delFilter 10000000  --realTimeLogging --logFile 16-fly-pg-${VERSION}.graphmap.log --maxCores 32 --outputFasta ./16-fly-pg-${VERSION}.gfa.fa.gz

cp cactus/src/cactus/cactus_progressive_config.xml ./config-split.xml
# hand-edit config-split.xml to have thresholds so as to be more permissive for smaller contig fragments.
# minQueryCoverages="0.25"
# minQueryCoverageThresholds=""

cactus-graphmap-split ./js ./softmask.fly.txt ./16-fly-pg-${VERSION}-minigraph.gfa.gz ./16-fly-pg-${VERSION}.paf --outDir ./chroms-16-fly-pg-${VERSION} --otherContig chrOther --refContigs $(for i in 2L 2R 3L 3R 4 X Y M; do echo chr$i; done) --reference dm6 --otherContig chrOther --realTimeLogging --logFile 16-fly-pg-${VERSION}.split.log --maxCores 32  --co
nfigFile ./config-split.xml

cactus-align-batch ./js ./chroms-16-fly-pg-${VERSION}/chromfile.txt align-16-fly-pg-${VERSION} --alignCores 8 --maxCores 32 --realTimeLogging --logFile 16-fly-pg-${VERSION}.align.log --alignOptions "--pangenome --pafInput --maxLen 10000 --reference dm6 --realTimeLogging  --outVG"

cactus-graphmap-join ./js --vg $(for j in 2L 2R 3L 3R 4 X Y M; do echo align-16-fly-pg-${VERSION}/chr$j.vg; done) --hal $(for j in 2L 2R 3L 3R 4 X Y M; do echo align-16-fly-pg-${VERSION}/chr$j.hal; done) --outDir join-16-fly-pg-${VERSION}-full --outName 16-fly-pg-${VERSION}-full --reference dm6 --gfaffix  --wlineSep "." --indexCores 31 --maxCores
 32 --realTimeLogging --logFile 16-fly-pg-${VERSION}.join-full.log

cactus-graphmap-join ./js --vg $(for j in 2L 2R 3L 3R 4 X Y M; do echo join-16-fly-pg-${VERSION}-full/clip-16-fly-pg-${VERSION}-full/chr$j.vg; done) --outDir join-16-fly-pg-${VERSION} --outName 16-fly-pg-${VERSION} --reference dm6  --wlineSep "." --clipLength 10000 --clipNonMinigraph --vcf --giraffe --preserveIDs --indexCores 32 --maxCores 32 --r
ealTimeLogging --logFile 16-fly-pg-${VERSION}.join.log

cactus-graphmap-join ./js --vg $(for j in 2L 2R 3L 3R 4 X Y M; do echo join-16-fly-pg-${VERSION}/clip-16-fly-pg-${VERSION}/chr$j.vg; done) --outDir join-16-fly-pg-${VERSION}-d2 --outName 16-fly-pg-${VERSION}-d2 --reference dm6  --wlineSep "." --vgClipOpts "-d 2 -m 1000" --giraffe --preserveIDs --indexCores 32 --maxCores 32 --realTimeLogging --log
File 16-fly-pg-${VERSION}.join-d2.log

```
