```
aws s3 cp s3://vg-k8s/users/hickey/fly-pangenome/16-fly-softmasked-fa.tar.gz .
tar zxf 16-fly-softmasked-fa.tar.gz

export VERSION=may26

cactus-minigraph ./js ./16-fly-pg-2022-05-26-seqfile.txt fly-pg-${VERSION}.gfa.gz --reference dm6 --maxCores 32 --mapCores 32 --realTimeLogging --logFile fly-pg-${VERSION}.minigraph.log

cactus-graphmap ./js ./16-fly-pg-2022-05-26-seqfile.txt ./fly-pg-${VERSION}.gfa.gz ./fly-pg-${VERSION}.paf --reference dm6 --delFilter 10000000  --realTimeLogging --logFile fly-pg-${VERSION}.graphmap.log --maxCores 32 --outputFasta ./fly-pg-${VERSION}.gfa.fa.gz

cp cactus/src/cactus/cactus_progressive_config.xml ./config-split.xml
# hand-edit config-split.xml to have thresholds so as to be more permissive for smaller contig fragments.
# minQueryCoverages="0.25"
# minQueryCoverageThresholds=""

cactus-graphmap-split ./js ./16-fly-pg-2022-05-26-seqfile.txt ./fly-pg-${VERSION}.gfa.gz ./fly-pg-${VERSION}.paf --outDir ./chroms-fly-pg-${VERSION} --otherContig chrOther --refContigs $(for i in 2L 2R 3L 3R 4 X Y M; do echo chr$i; done) --reference dm6 --otherContig chrOther --realTimeLogging --logFile fly-pg-${VERSION}.split.log --maxCores 32  --co
nfigFile ./config-split.xml

cactus-align-batch ./js ./chroms-fly-pg-${VERSION}/chromfile.txt align-fly-pg-${VERSION} --alignCores 8 --maxCores 32 --realTimeLogging --logFile fly-pg-${VERSION}.align.log --alignOptions "--pangenome --pafInput --maxLen 10000 --reference dm6 --realTimeLogging  --outVG"

cactus-graphmap-join ./js --vg $(for j in 2L 2R 3L 3R 4 X Y M; do echo align-fly-pg-${VERSION}/chr$j.vg; done) --hal $(for j in 2L 2R 3L 3R 4 X Y M; do echo align-fly-pg-${VERSION}/chr$j.hal; done) --outDir join-fly-pg-${VERSION}-full --outName fly-pg-${VERSION}-full --reference dm6 --gfaffix  --wlineSep "." --indexCores 31 --maxCores
 32 --realTimeLogging --logFile fly-pg-${VERSION}.join-full.log

cactus-graphmap-join ./js --vg $(for j in 2L 2R 3L 3R 4 X Y M; do echo join-fly-pg-${VERSION}-full/clip-fly-pg-${VERSION}-full/chr$j.vg; done) --outDir join-fly-pg-${VERSION} --outName fly-pg-${VERSION} --reference dm6  --wlineSep "." --clipLength 10000 --clipNonMinigraph --vcf --giraffe --preserveIDs --indexCores 32 --maxCores 32 --r
ealTimeLogging --logFile fly-pg-${VERSION}.join.log

cactus-graphmap-join ./js --vg $(for j in 2L 2R 3L 3R 4 X Y M; do echo join-fly-pg-${VERSION}/clip-fly-pg-${VERSION}/chr$j.vg; done) --outDir join-fly-pg-${VERSION}-d2 --outName fly-pg-${VERSION}-d2 --reference dm6  --wlineSep "." --vgClipOpts "-d 2 -m 1000" --giraffe --preserveIDs --indexCores 32 --maxCores 32 --realTimeLogging --log
File fly-pg-${VERSION}.join-d2.log

```
