library(dplyr) ## install with: install.packages('dplyr')
library(ggplot2) ## install with: install.packages('ggplot2')

## read dm6-bwa mapping stats
map.d = list.files('dm6_mappings/')
map.d = lapply(map.d, function(ff){
  df = read.table(paste0('dm6_mappings/', ff), as.is=TRUE)
  colnames(df) = c('n', 'mapq', 'perfect')
  df$sample = gsub('dm6.(.*).bwa.mapstats.txt', '\\1', ff)
  df
}) %>% bind_rows

## read pangenome-giraffe mapping stats
map.g = list.files('fly-pg-may26-d2_mappings/')
map.g = lapply(map.g, function(ff){
  df = read.table(paste0('fly-pg-may26-d2_mappings/', ff), as.is=TRUE)
  colnames(df) = c('n', 'mapq', 'perfect')
  df$sample = gsub('fly-pg-may26-d2.(.*).giraffe.mapstats.txt', '\\1', ff)
  df
}) %>% bind_rows

## combine
map.df = rbind(map.d %>% mutate(method='dm6_bwa'),
               map.g %>% mutate(method='cactus_giraffe')
               ) %>%
  mutate(perfect=as.logical(perfect))

## subset to samples analyzed by both mappers
map.df = map.df %>% group_by(sample) %>% filter(length(unique(method))>1)

map.df = map.df %>% mutate(method=factor(method, levels=c('dm6_bwa', 'cactus_giraffe'),
                                         labels=c('dm6-BWA', 'Cactus-Giraffe')))

gp.df = map.df %>% 
  group_by(sample, method, perfect) %>%
  summarize(n=sum(n)) %>% 
  group_by(sample, method) %>% mutate(prop=n/sum(n)) %>%
  filter(perfect)
gp.df %>% group_by(method) %>% summarize(prop=mean(prop))
ggp.perfect = ggplot(gp.df, aes(x=method, y=prop)) +
  geom_line(aes(group=sample), alpha=.1) + 
  geom_boxplot(alpha=.4, fill='lightblue') + 
  theme_bw() +
  coord_flip() + 
  ylab('proportion of reads aligned perfectly')
ggp.perfect

gp.df = map.df %>% 
  mutate(mapping.status=ifelse(mapq>0, 'high', 'low')) %>%
  group_by(sample, method, mapping.status) %>%
  summarize(n=sum(n)) %>% 
  group_by(sample, method) %>% mutate(prop=n/sum(n)) %>%
  filter(mapping.status=='high')
gp.df %>% group_by(method) %>% summarize(prop=mean(prop))
ggp.mapq = ggplot(gp.df, aes(x=method, y=prop)) +
  geom_line(aes(group=sample), alpha=.1) + 
  geom_boxplot(alpha=.3, fill='indianred2') + 
  theme_bw() + coord_flip() + 
  ylab('proportion of reads with mapq>0')
ggp.mapq

pdf('fly-mapping-stats.pdf', 8, 2)
ggp.perfect
ggp.mapq
dev.off()
