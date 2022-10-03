library(dplyr)
library(ggplot2)
library(ggrepel)

args = commandArgs(TRUE)
## Inputs: summaryc csv path/file names

df = lapply(args, function(f.path){
  fn = basename(f.path)
  read.csv(f.path, as.is=TRUE) %>%
    mutate(file=fn,
           truthset=gsub('(.*)\\..*\\..*\\..*\\..*\\.roc\\.all\\.csv\\.gz', '\\1', file),
           sample=gsub('.*\\.(.*)\\..*\\..*\\..*\\.roc\\.all\\.csv\\.gz', '\\1', file),
           reads=gsub('.*\\..*\\.(.*)\\..*\\..*\\.roc\\.all\\.csv\\.gz', '\\1', file),
           method=gsub('.*\\..*\\..*\\.(.*)\\..*\\.roc\\.all\\.csv\\.gz', '\\1', file),
           region=gsub('.*\\..*\\..*\\..*\\.(.*)\\.roc\\.all\\.csv\\.gz', '\\1', file)) %>%
    select(-file)
}) %>% bind_rows

df = df %>% filter(Subtype=='*', Filter=='ALL', Subset=="*", QQ!='*') %>%
  mutate(QQ=as.numeric(QQ), precision=METRIC.Precision, recall=METRIC.Recall, f1=METRIC.F1_Score) %>%
  select(truthset, sample, reads, method, region, Type, QQ, precision, recall, f1, method, TRUTH.TOTAL, TRUTH.TP, TRUTH.FN, QUERY.FP, QUERY.TP, QUERY.TOTAL)

write.table(df, file='eval-roc-summary.tsv', sep='\t', row.names=FALSE, quote=FALSE)

pdf('eval-roc-summary.pdf', 9, 6)
df.best = df %>% group_by(Type, method, reads, sample, region, truthset) %>% arrange(desc(f1)) %>% do(head(., 1)) %>%
  mutate(label=paste0('F1:', f1))

tmp = lapply(unique(df$regions), function(reg){
  tmp = lapply(unique(df$reads), function(read){
  
    ggp = df %>% filter(reads==read, region==reg) %>% 
      arrange(method, Type, QQ) %>%
      ggplot(aes(x=precision, y=recall, color=method)) +
      geom_path(size=1.5, alpha=.8) +
      theme_bw() +
      facet_wrap(paste(sample, truthset)~Type, scales='free') +
      theme(legend.position='bottom') +
      scale_color_brewer(palette='Set1') +
      geom_label_repel(aes(label=label), data=df.best, show.legend=FALSE) +
      ggtitle(reg) + labs(caption=read)
    print(ggp)
    
  })  
})
dev.off()
