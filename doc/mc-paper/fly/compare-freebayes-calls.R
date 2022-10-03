library(dplyr) ## install with: install.packages('dplyr')
library(ggplot2) ## install with: install.packages('ggplot2')

## read the comparison files for all samples
files = list.files('freebayes_stats', 'stats.tsv.gz')
df = lapply(files, function(ff){
  df = read.table(paste0('freebayes_stats/', ff), as.is=TRUE, sep='\t', header=TRUE)
  df %>% mutate(sample=gsub("(.*).isec.stats.tsv.gz", '\\1', ff))
}) %>% bind_rows

## using the QUAL field as quality threshod
df.qual = df %>% filter(qual.type=='QUAL') %>% arrange(min.qual) %>%
  mutate(qual.c=factor(qual.c, levels=unique(qual.c))) %>%
  filter(call!='hom ref') %>% group_by(set, sample) %>% mutate(prop=n/sum(n))

## using the GQ field as quality threshod
df.gq = df %>% filter(qual.type=='GQ') %>% arrange(min.qual) %>%
  mutate(qual.c=factor(qual.c, levels=unique(qual.c))) %>%
  filter(call!='hom ref') %>% group_by(set, sample) %>% mutate(prop=n/sum(n))

## figures
pdf('freebayes-comp.pdf', 8, 4)

ggp.q = df.qual %>% filter(min.qual>.1) %>%
  group_by(sample, set) %>% summarize(n=sum(n)) %>%
  mutate(set.c=ifelse(grepl('only', set), 'unique to a method', 'called by both methods'),
         method=ifelse(grepl('bwa', set), 'BWA on dm6', 'vg giraffe on\ncactus pangenome')) %>%
  ggplot(aes(x=n/1000, fill=method)) +
  geom_histogram(position='dodge') +
  theme_bw() +
  facet_grid(.~set.c, scales='free') +
  scale_fill_brewer(palette='Set1', name='mapping approach') +
  xlab('number of variant sites (K)') +
  ylab('number of samples')
ggp.q

ggp.q.c = lapply(10^(-7:2), function(th){
  df.qual %>% filter(min.qual>th) %>%
    group_by(sample, set) %>% summarize(n=sum(n)) %>%
    mutate(min.qual.th=th)
}) %>% bind_rows %>% 
  mutate(set.c=ifelse(grepl('only', set), 'unique to a method', 'called by both methods'),
         method=ifelse(grepl('bwa', set), 'BWA on dm6', 'vg giraffe on\ncactus pangenome')) %>%
  ggplot(aes(x=min.qual.th, y=n/1000, fill=method)) +
  geom_boxplot(aes(group=paste(method, min.qual.th))) +
  theme_bw() +
  ylab('number of variant sites (K)') +
  xlab('minimum quality threshold (QUAL)') +
  facet_wrap(~set.c, scales='free', nrow=1) +
  scale_fill_brewer(palette='Set1', name='mapping approach') +
  scale_x_log10(breaks=10^(-7:2)) +
  theme(legend.position='bottom', axis.text.x=element_text(angle=30, hjust=1))
ggp.q.c

ggp.gq = df.gq %>% filter(min.qual>10) %>%
  group_by(sample, set) %>% summarize(n=sum(n)) %>%
  mutate(set.c=ifelse(grepl('only', set), 'unique to a method', 'called by both methods'),
         method=ifelse(grepl('bwa', set), 'BWA on dm6', 'vg giraffe on\ncactus pangenome')) %>%
  ggplot(aes(x=n/1000, fill=method)) +
  geom_histogram(position='dodge') +
  theme_bw() +
  facet_grid(.~set.c, scales='free') +
  scale_fill_brewer(palette='Set1', name='mapping approach') +
  xlab('number of variant sites (K)') +
  ylab('number of samples')
ggp.gq

ggp.gq.c = lapply(seq(0,160,10), function(th){
  df.gq %>% filter(min.qual>th) %>%
    group_by(sample, set) %>% summarize(n=sum(n)) %>%
    mutate(min.qual.th=th)
}) %>% bind_rows %>% 
  mutate(set.c=ifelse(grepl('only', set), 'unique to a method', 'called by both methods'),
         method=ifelse(grepl('bwa', set), 'BWA on dm6', 'vg giraffe on\ncactus pangenome')) %>%
  ggplot(aes(x=min.qual.th, y=n/1000, fill=method)) +
  geom_boxplot(aes(group=paste(method, min.qual.th))) +
  theme_bw() +
  facet_wrap(~set.c, scales='free') +
  ylab('number of variant sites (K)') +
  xlab('minimum genotype quality (GQ)') +
  scale_fill_brewer(palette='Set1', name='mapping approach') + 
  theme(legend.position='bottom')
ggp.gq.c

dev.off()

## numbers
df.qual %>% filter(min.qual>.1) %>%
  mutate(method=ifelse(grepl('bwa', set), 'bwa', 'giraffe')) %>%
  group_by(method, sample) %>% summarize(n=sum(n)) %>%
  group_by(method) %>% summarize(n=mean(n))

df.gq %>% filter(min.qual>10) %>%
  mutate(method=ifelse(grepl('bwa', set), 'bwa', 'giraffe')) %>%
  group_by(method, sample) %>% summarize(n=sum(n)) %>%
  group_by(method) %>% summarize(n=mean(n))
