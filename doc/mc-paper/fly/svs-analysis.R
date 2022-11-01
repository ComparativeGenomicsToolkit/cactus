library(dplyr) ## install with: install.packages('dplyr')
library(ggplot2) ## install with: install.packages('ggplot2')
library(sveval) ## install with: BiocManager::install('jmonlong/sveval')
## install BiocManager with: install.packages('BiocManager')

## how to split SVs in to size classes in the graphs
size.breaks = c(-1, 0, 10, 49, 100, 500, 1000, 5000, 1e4, 1e5, 1e6,Inf)
size.labs = c('0', '1-10', '11-49', '50-100', '100-500',
              '500-1K', '1K-5K', '5K-10K', '10K-100K', '100K-1M', '>1M')

## SVs in the pangenome
svs = readRDS('16-fruitfly-mc-2022-05-26.svs.site.rol90.insd100.rds')

n.samps = length(unique(svs$sample))

svs.df = svs %>% as.data.frame %>%
  group_by(svsite, type, clique, seqnames, start, end, size) %>%
  mutate(n=length(unique(sample))) %>%
  group_by(svsite, type, clique, seqnames) %>%
  arrange(desc(n)) %>% 
  summarize(start=min(start), end=max(end),
    size.min=min(size), size.max=max(size),
    size=size[1], alleles=length(unique(paste(start, end, size, alt, ref))),
    ac=length(unique(sample)), af=ac/n.samps,
    .groups='drop')

## calls across cohort of 100 samples

calls = readRDS('16-fruitfly-mc-2022-05-26-d2.100samples.decomposed.svs.site.rol90.insd100.rds')

calls.df = calls %>% as.data.frame %>%
  filter(size>=40) %>% 
  group_by(svsite, type, clique, seqnames) %>%
  arrange(desc(af)) %>% 
  summarize(start=min(start), end=max(end),
    size.min=min(size), size.max=max(size),
    size=size[1], alleles=n(),
    ac=sum(ac), af=sum(af),
    .groups='drop') %>%
  mutate(af=ifelse(af>1, 1, af))

## add RepeatMasker annotation
calls.rm = readRDS('16-fruitfly-mc-2022-05-26-d2.100samples.decomposed.svs.site.rol90.insd100.rmsk.rds')
calls.rm = calls.rm %>% as.data.frame %>% select(svsite, rmsk.classfam, rmsk.name, rmsk.cov) %>% unique
calls.df = merge(calls.df, calls.rm, all.x=TRUE)

## overlap SVs in pangenome with calls
ol.o = svOverlap(makeGRangesFromDataFrame(calls.df, keep.extra.columns=TRUE),
                 makeGRangesFromDataFrame(svs.df, keep.extra.columns=TRUE),
                 min.ol=.9,
                 max.ins.dist=100)
ol.o.svs = ol.o %>% as.data.frame %>% mutate(calls.af=calls.df$af[queryHits]) %>%
  group_by(subjectHits) %>% summarize(calls.af=mean(calls.af))

svs.df$calls.af = 0
svs.df$calls.af[ol.o.svs$subjectHits] = ol.o.svs$calls.af

## figures

ggp.freq = ggplot(svs.df, aes(x=ac, fill=type)) +
  geom_bar(position=position_dodge(preserve='single'), color='black', size=.4) +
  theme_bw() + 
  ylab('number of SV sites') +
  xlab('alternate allele count in pangenome') +
  scale_fill_brewer(name='SV type', palette='Set2') + 
  scale_x_continuous(breaks=0:15) +
  scale_y_log10(breaks=1+10^(0:10), labels=10^(0:10)) +
  theme(legend.position=c(.99,.99), legend.justification=c(1,1)
        )
ggp.freq

ggp.size = svs.df %>% mutate(size=cut(size, breaks=size.breaks, labels=size.labs),
                             called=ifelse(calls.af>0, 'called in 100\nshort-reads samples', 'only in pangenome'),
                             called=factor(called, levels=rev(sort(unique(called))))) %>%
  group_by(size, type, called) %>% summarize(n=n()) %>%
  ggplot(aes(x=size, y=n, fill=called)) +
  geom_bar(stat='identity', color='black', size=.4) +
  xlab('size') +
  ylab('number of SV sites') +
  theme_bw() +
  facet_wrap(~type, scales='free_y', nrow=1) +
  xlab('SV size (bp)') + 
  scale_fill_brewer(palette='Greys') +
  theme(legend.title=element_blank(),
        axis.text.x=element_text(angle=50, hjust=1))
ggp.size

af.breaks = c(-1,0,.1,.3,.6,1)
af.labels = c(0, '0-0.1', '0.1-0.3', '0.3-0.6', '0.6-1')
ggp.call.freq = svs.df %>% 
  mutate(ac=cut(ac, c(0,1,2,3,Inf), labels=c(1:3, '4+')),
         calls.af=cut(calls.af, af.breaks, af.labels)) %>% 
  ggplot(aes(x=ac, fill=calls.af)) +
  geom_bar(color='black', size=.4) +
  ylab('number of SV sites') +
  facet_wrap(.~ac, scales='free', nrow=1) + 
  theme_bw() +
  xlab('alternate allele count in pangenome') +
  scale_fill_brewer(name='allele frequency\nin the cohort\nof 100 samples') +
  theme(strip.text=element_blank())
ggp.call.freq

ggp.freq.size = calls.df %>%
  mutate(af=cut(af, af.breaks, af.labels)) %>%
  ggplot(aes(x=size, fill=af)) +
  geom_histogram(bins=90, color='black') +
  scale_fill_brewer(name='allele frequency') +
  ylab('number of SV sites') +
  xlab('SV size (bp)') + 
  theme_bw() +
  scale_x_log10() +
  theme(legend.position=c(.99,.99), legend.justification=c(1,1))
ggp.freq.size

## function to keep the top N most common values and group the rest in 'other'
## used to highlight the most common repeat classes
factop <- function(x, top.n=5){
  xt = head(sort(table(x), decreasing=TRUE), top.n)
  factor(ifelse(x %in% names(xt), x, 'other'), levels=c(names(xt), 'other'))
}

ggp.rep.size = calls.df %>%
  mutate(repeat.type=ifelse(!is.na(rmsk.cov) & rmsk.cov>.8, rmsk.classfam, 'none'),
         repeat.type=factop(repeat.type, 7)) %>%
  ggplot(aes(x=size, fill=repeat.type)) +
  geom_histogram(bins=90, color='black') +
  ylab('number of SV sites') +
  xlab('SV size (bp)') + 
  theme_bw() +
  scale_fill_brewer(palette='Set2', name='repeat type') + 
  scale_x_log10() +
  theme(legend.position=c(.99,.99), legend.justification=c(1,1))
ggp.rep.size

## make a PDF with figures
pdf('16-fruitfly-mc-2022-05-26.svs.pdf', 9, 4)
ggp.freq
ggp.size
ggp.call.freq
ggp.freq.size
ggp.rep.size
ggp.repn.size
dev.off()

## quick look at inversions
svs.df %>% filter(type=='INV') %>% select(seqnames, start, end, size) %>% arrange(desc(size))

