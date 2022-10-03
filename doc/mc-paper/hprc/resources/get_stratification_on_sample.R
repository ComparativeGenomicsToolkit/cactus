args = commandArgs(TRUE)
## 1. tsv made with: bcftools query -f '%CHROM\t%POS\t%END[\t%BD\t%BVT]\n' happy-annotated.vcf.gz | awk '{if(($4!="." && $4 != "N") || ($6!="." && $6 != "N")){print $0}}' | gzip
## 2. output RDS file

library(dplyr)
library(GenomicRanges)
library(tidyr)

## summary counts
sumEval <-function(vars){
  ss.df = mcols(vars) %>% as.data.frame %>% group_by(eval, type) %>% summarize(n=n(), .groups='drop') %>%
    pivot_wider(names_from=eval, values_from=n, values_fill=0)
  if(is.null(ss.df$FN)) ss.df$FN = 0
  if(is.null(ss.df$TP)) ss.df$TP = 0
  if(is.null(ss.df$FP)) ss.df$FP = 0
  ss.df %>% 
    mutate(called=TP+FP,
           precision=TP/(TP + FP), precision=round(precision, 6),
           recall=TP/(TP + FN), recall=round(recall, 6),
           recall=ifelse(is.nan(recall), 0, recall),
           F1 = 2 * precision * recall/(precision + recall), F1 = round(F1, 6),
           F1 = ifelse(precision == 0 & recall == 0, 0, F1)) %>%
    select(type, everything())
}
## equivalent GRanges object for overlap later
extractVariants <- function(evm){
  evm.gr = makeGRangesFromDataFrame(evm)
  tp.indel = evm.gr[which(evm$giab=='TP' & evm$giab.type=='INDEL')]
  tp.indel$type = 'INDEL'
  tp.indel$eval = 'TP'
  fn.indel = evm.gr[which(evm$giab=='FN' & evm$giab.type=='INDEL')]
  fn.indel$type = 'INDEL'
  fn.indel$eval = 'FN'
  fp.indel = evm.gr[which(evm$grch38=='FP' & evm$grch38.type=='INDEL')]
  fp.indel$type = 'INDEL'
  fp.indel$eval = 'FP'
  tp.snp = evm.gr[which(evm$giab=='TP' & evm$giab.type=='SNP')]
  tp.snp$type = 'SNP'
  tp.snp$eval = 'TP'
  fn.snp = evm.gr[which(evm$giab=='FN' & evm$giab.type=='SNP')]
  fn.snp$type = 'SNP'
  fn.snp$eval = 'FN'
  fp.snp = evm.gr[which(evm$grch38=='FP' & evm$grch38.type=='SNP')]
  fp.snp$type = 'SNP'
  fp.snp$eval = 'FP'
  return(c(tp.indel, fn.indel, fp.indel, tp.snp, fn.snp, fp.snp))
}

## cores to use
NB.CORES = 2

## read happy results
eval = read.table(args[1], as.is=TRUE, sep='\t', header=FALSE)
colnames(eval) = c('chr', 'start', 'end', 'giab', 'giab.type', 'grch38', 'grch38.type')

vars = extractVariants(eval)

giab.regs = readRDS('regions.all.giab.stratifications.RDS')

eval.giab.regs = mclapply(names(giab.regs), function(regn){
  ol = overlapsAny(vars, giab.regs[[regn]])
  sumEval(vars[which(ol),]) %>% mutate(region=regn)
}, mc.cores=NB.CORES) %>% bind_rows

saveRDS(eval.giab.regs, file=args[2])
