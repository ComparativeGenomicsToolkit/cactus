This is the same input as `5-cow-pg-2022-09-22` but run with the new `cactus-pangenome` interface to do everything in one command. 

Cactus commit:  b5982ed1fa54aaf4c22087ab4132098230b97fe3

Note: If I could go back, I'd probably add a `--gbz` to get a GBZ output of the clipped graph (can be obtained here using `vg gbwt` on the `gfa` output)

```
cactus-pangenome aws:us-west-2:glennhickey-jobstore-cow-pg ./5-cow-pg-2023-03-31-seqfile.txt --outDir s3://vg-k8s/users/hickey/5-cow-pg-2023-03-31 --outName 5-cow-pg-2023-03-31 --indexCores 31 --mapCores 8 --indexCores 31 --consCores 16 --logFile 5-cow-pg-2023-03-31.log  --giraffe --vcf --gfa  --otherContig chrOther --refContigs $(for i in `seq 29`; do echo chr$i; done ; echo "chrX chrM") --reference bosTau9 --batchSystem mesos --provisioner aws --defaultPreemptable --nodeType r5.8xlarge:1.5 --nodeStorage 650 --maxNodes 25 --betaInertia 0 --targetTime 1
```