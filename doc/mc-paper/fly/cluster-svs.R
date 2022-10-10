library(sveval) ## install with: BiocManager::install('jmonlong/sveval')
library(GenomicRanges) ## install with: BiocManager::install('GenomicRanges')
## install BiocManager with: install.packages('BiocManager')

## simple repeat annotation
if(!file.exists('dm6.simpleRepeat.txt.gz')){
  download.file('https://hgdownload.soe.ucsc.edu/goldenPath/dm6/database/simpleRepeat.txt.gz', 'dm6.simpleRepeat.txt.gz')
}
sr = read.table('dm6.simpleRepeat.txt.gz', as.is=TRUE, sep='\t')
sr = GRanges(paste0('dm6.', sr$V2), IRanges(sr$V3, sr$V4))
sr = reduce(sr)

## SVs in the pangenome
svs = lapply(list.files('svs', '*.rds'), function(fn){
  gr = readRDS(paste0('svs/', fn))
  gr$sample = gsub('.+\\.(.*)\\.decomposed.svs.rds', '\\1', fn)
  gr
})

svs = clusterSVs(svs, range.seq.comp=TRUE, ins.seq.comp=TRUE, nb.cores=16, batch.maxsize=500, simprep=sr, min.rol=.9, max.ins.dist=100)

saveRDS(svs, file='fly-pg-may26.svs.site.rol90.insd100.rds')

## calls
svs = readRDS('fly-pg-may26-d2.100samples.decomposed.svs.rds')
svs = subset(svs, ac>0)

svs = clusterSVs(svs, range.seq.comp=TRUE, ins.seq.comp=TRUE, nb.cores=16, batch.maxsize=500, simprep=sr, min.rol=.9, max.ins.dist=100)

saveRDS(svs, file='fly-pg-may26-d2.100samples.decomposed.svs.site.rol90.insd100.rds')
