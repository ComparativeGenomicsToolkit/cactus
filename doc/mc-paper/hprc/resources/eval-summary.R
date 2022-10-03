library(dplyr)
library(ggplot2)

args = commandArgs(TRUE)
## Inputs: summaryc csv path/file names

df = lapply(args, function(f.path){
  fn = basename(f.path)
  read.csv(f.path, as.is=TRUE) %>%
    mutate(file=fn,
           truthset=gsub('(.*)\\..*\\..*\\..*\\..*\\.summary\\.csv', '\\1', file),
           sample=gsub('.*\\.(.*)\\..*\\..*\\..*\\.summary\\.csv', '\\1', file),
           reads=gsub('.*\\..*\\.(.*)\\..*\\..*\\.summary\\.csv', '\\1', file),
           method=gsub('.*\\..*\\..*\\.(.*)\\..*\\.summary\\.csv', '\\1', file),
           region=gsub('.*\\..*\\..*\\..*\\.(.*)\\.summary\\.csv', '\\1', file)) %>%
    select(-file)
}) %>% bind_rows

write.table(df, file='eval-summary.tsv', sep='\t', row.names=FALSE, quote=FALSE)

pdf('eval-summary.pdf', 9, 6)
for(metric in c("METRIC.F1_Score", "METRIC.Recall", "METRIC.Precision", 'TRUTH.TP', 'TRUTH.FN', 'QUERY.FP')){
  for(reg in unique(df$region)){
    
    cat('\n\n## ', metric, ' - ', reg, '\n\n')
  
  ggp = df %>% filter(Filter=='PASS', region==reg) %>% 
    ggplot(aes_string(x='sample', color='method', y=metric)) +
    geom_point(stat='identity', position=position_dodge(.5)) + 
    facet_grid(Type~truthset, space='free', scales='free') +
    theme_bw()
  print(ggp)
  }  
}
def.off()
