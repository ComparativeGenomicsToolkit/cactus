library(GenomicRanges) ## install with: BiocManager::install('GenomicRanges')
## install BiocManager with: install.packages('BiocManager')
library(knitr) ## install with: install.packages('knitr')

rds.f = list.files('happy_runs/', 'stratification.rds')

## read the evaluation in each region set for each sample and method
eval.giab.regs = lapply(rds.f, function(fn){
  df = readRDS(paste0('happy_runs/', fn))
  df %>% mutate(sample=gsub('.+\\.(.+)\\..+\\..+\\.nofalsedup_in_chm13.stratification.rds', '\\1', fn),
                method=gsub('.+\\..+\\..+\\.(.+)\\.nofalsedup_in_chm13.stratification.rds', '\\1', fn))
}) %>% bind_rows

## compute numbers for all samples
eval.giab.regs = eval.giab.regs %>%
  mutate(method=factor(method,
                       levels=c("giraffedv", "giraffedv_chm13"),
                       labels=c('HPRC_GRCh38', 'HPRC_CHM13'))) %>% 
  group_by(method, region) %>% 
  summarize(TP=sum(TP, na.rm=TRUE),
            FP=sum(FP, na.rm=TRUE),
            FN=sum(FN, na.rm=TRUE)) %>%
  mutate(precision=TP/(TP+FP), recall=TP/(TP + FN),
         F1=2*precision * recall/(precision + recall))

## amount of sequence for each region set
giab.regs = readRDS('regions.all.giab.stratifications.RDS')
regs.mb = tibble(region=names(giab.regs),
                 mb=unlist(lapply(giab.regs, function(gr) sum(width(gr)/1e6))))

## those regions were excluded from the evaluation
regions.torm = c("false_duplications_correct_copy",
                 "false_duplications_incorrect_copy",
                 "collapsed_duplication_FP_regions",
                 "CMRGv1.00_falselyduplicatedgenes",
                 "population_CNV_FP_regions")

pdf('hprc-smallvariant-giab-hg00125-stratifications.pdf', 8, 8)
eval.giab.regs %>% filter(!(region %in% regions.torm)) %>%
  dplyr::select(method, region, F1) %>% 
  pivot_wider(names_from=method, values_from=F1) %>%
  mutate(f1.diff=HPRC_CHM13-HPRC_GRCh38) %>%
  merge(regs.mb) %>%
  ungroup %>%
    mutate(label=ifelse(rank(-abs(f1.diff))<15, region, '')) %>% 
    ggplot(aes(y=f1.diff, x=mb*1e6)) +
    geom_hline(yintercept=0, linetype=2, alpha=.7) + 
    geom_point() +
    geom_label_repel(aes(label=label), show.legend=FALSE) + 
    theme_bw() +
    scale_x_log10() +
    xlab('total size (bp)') +
    ylab('F1 score difference (CHM13-GRCh38)')
dev.off()

