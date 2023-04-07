library(dplyr) ## install with: install.packages('dplyr')
library(ggplot2) ## install with: install.packages('ggplot2')
library(RColorBrewer)
library(colorspace)

## color palette
pal = brewer.pal(4, 'Set1')[c(4,1,1)]
pal[2] = darken(pal[2], .3)
pal[3] = lighten(pal[3], .2)

## snarls
snarls.df = rbind(
  read.table('hprc-v1.0-mc-grch38-maxdel.10mb.dist-stats.tsv.gz',
             as.is=TRUE, header=TRUE) %>% mutate(pangenome='minigraph-cactus'),
  read.table('hprc-v1.0-mc-chm13-maxdel.10mb.dist-stats.tsv.gz',
             as.is=TRUE, header=TRUE) %>% mutate(pangenome='minigraph-cactus-chm13'),
  read.table('hprc-v1.0-minigraph-grch38.dist-stats.tsv.gz',
             as.is=TRUE, header=TRUE) %>% mutate(pangenome='minigraph')
)

snarls.df = snarls.df %>% mutate(diff_length=max_length-min_length)

## how to split snarls in to size classes in the graphs
size.breaks = c(-1, 0, 10, 49, 100, 500, 1000, 5000, 1e4, 1e5, Inf)
size.labs = c('0', '1-10', '11-49', '50-100', '100-500',
              '500-1K', '1K-5K', '5K-10K', '10K-100K', '>100K')
y.bks = 10^(1:10)

snarls.df %>%
  mutate(diff_length=cut(diff_length, breaks=size.breaks, labels=size.labs),
         pangenome=factor(pangenome, levels=c('minigraph', 'minigraph-cactus',
                                              'minigraph-cactus-chm13'),
                          labels=c('Minigraph', 'GRCh38_Minigraph-Cactus',
                                   'CHM13_Minigraph-Cactus'))) %>%
  group_by(pangenome, diff_length) %>% summarize(n=sum(n)) %>%
  ggplot(aes(x=diff_length, y=n+1, fill=pangenome)) +
  geom_bar(stat='identity', position='dodge') +
  geom_text(aes(label=paste0(format(n, big.mark=','), '  ')), angle=90,
            position=position_dodge(.9), hjust=1, color='white', size=2) + 
  xlab('difference between shortest and largest paths through snarl (bp)') +
  ylab('number of snarls') +
  theme_bw() +
  scale_fill_manual(values=pal) + 
  theme(legend.position=c(.99,.99), legend.justification=c(1,1),
        legend.title=element_blank()) + 
  scale_y_log10(breaks=y.bks, labels=format(y.bks, big.mark=',', scientific=FALSE))
