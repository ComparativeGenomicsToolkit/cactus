library(dplyr)
library(GenomicRanges)
library(ggplot2)

if(!file.exists('v3.0-GRCh38-stratifications-all-except-genome-specific-stratifications.tsv')){
  download.file('https://raw.githubusercontent.com/genome-in-a-bottle/genome-stratifications/master/GRCh38/v3.0-GRCh38-stratifications-all-except-genome-specific-stratifications.tsv', 'v3.0-GRCh38-stratifications-all-except-genome-specific-stratifications.tsv')
}
reg.tsv = read.table('v3.0-GRCh38-stratifications-all-except-genome-specific-stratifications.tsv', as.is=TRUE, sep='\t')
colnames(reg.tsv) = c('name', 'path')

regions.l = lapply(1:nrow(reg.tsv), function(ii){
  download.file(paste0('https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/genome-stratifications/v3.0/GRCh38/', reg.tsv$path[ii]), 'temp.bed.gz')
  reg = read.table('temp.bed.gz', as.is=TRUE, sep='\t')
  reg = GRanges(reg[,1], IRanges(reg[,2], reg[,3]))
  unlink('temp.bed.gz')
  gc()
  return(reg)
})
names(regions.l) = reg.tsv$name
saveRDS(regions.l, file='regions.all.giab.stratifications.RDS')
