library(dplyr) ## install with: install.packages('dplyr')
library(ggplot2) ## install with: install.packages('ggplot2')
library(tidyr) ## install with: install.packages('tidyr')
library(RColorBrewer)
library(knitr)  ## install with: install.packages('knitr')

map.f = list.files('mapstats/')
map.all = lapply(map.f, function(ff){
  df = read.table(paste0('mapstats/', ff), as.is=TRUE)
  colnames(df) = c('n', 'mapq', 'perfect')
  df$method = ifelse(grepl('giraffe', ff), 'cactus_giraffe', 'grch38_bwamem')
  df$method = ifelse(grepl('chm13', ff), 'cactus_giraffe_chm13', df$method)
  df$method = ifelse(grepl('noclip', ff), 'cactus_giraffe_noclip', df$method)
  df$sample = gsub('([^\\.]*)\\..*', '\\1', ff)
  df
}) %>% bind_rows

## all methods
meth.levels = c('grch38_bwamem', 'cactus_giraffe_noclip', 'cactus_giraffe',
                'cactus_giraffe_chm13')
meth.labels = c('GRCh38-BWA', 'Cactus_GRCh38_noclipping-Giraffe', 'Cactus_GRCh38-Giraffe',
                'Cactus_CHM13-Giraffe')

## the graphs can include the full graph (no "clipping") or not based on this value
show.clip = FALSE
## the main figure doesn't include stats for this graph
## to make the supp figure that does, use `show.clip = TRUE`

cbbPalette = c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
col.pal = brewer.pal(4, 'Set2')[c(1,3,2)]
if(!show.clip){
  meth.levels = c('grch38_bwamem', 'cactus_giraffe',
                  'cactus_giraffe_chm13')
  meth.labels = c('GRCh38-BWA', 'HPRC_GRCh38-Giraffe', 'HPRC_CHM13-Giraffe')
  out.pdf = ''
} else {
  ## subset for manuscript
  meth.levels = c('grch38_bwamem', 'cactus_giraffe_noclip', 'cactus_giraffe',
                  'cactus_giraffe_chm13')
  meth.labels = c('GRCh38-BWA', 'HPRC_GRCh38_full-Giraffe', 'HPRC_GRCh38-clipped-Giraffe', 'HPRC_CHM13-clipped-Giraffe')
  out.pdf = '-noclip'
  col.pal = c(col.pal[1], brewer.pal(4, 'Paired')[2:1], col.pal[3])
}

map.df = map.all %>% 
  mutate(perfect=as.logical(perfect),
         method=factor(method, levels=meth.levels, labels=meth.labels)) %>%
  filter(!is.na(method))

gp.df = map.df %>% 
  mutate(mapping.status=ifelse(mapq==-1, 'unmapped', 'mapped')) %>%
  group_by(sample, method, mapping.status) %>%
  summarize(n=sum(n)) %>% 
  group_by(sample, method) %>% mutate(prop=n/sum(n)) %>%
  filter(mapping.status=='unmapped')
gp.df %>% group_by(method) %>% summarize(prop=mean(prop))
ggp.map = ggplot(gp.df, aes(x=sample, y=prop, colour=method)) +
  geom_point(position=position_dodge(.3)) + 
  theme_bw() +
  scale_color_manual(values=col.pal) + 
  coord_flip() + 
  ylab('proportion of unmapped reads')
ggp.map

gp.df = map.df %>% 
  group_by(sample, method, perfect) %>%
  summarize(n=sum(n)) %>% 
  group_by(sample, method) %>% mutate(prop=n/sum(n)) %>%
  filter(perfect)
gp.df %>% group_by(method) %>% summarize(prop=mean(prop))
ggp.perfect = ggplot(gp.df, aes(x=sample, y=prop, fill=method)) +
  geom_bar(stat='identity', position=position_dodge(.8), width=.7) + 
  scale_fill_manual(values=col.pal) + 
  theme_bw() + coord_flip() + 
  ylab('proportion of reads aligned perfectly') +
  ylim(0,1)
ggp.perfect

gp.df = map.df %>% 
  mutate(mapping.status=ifelse(mapq>0, 'high', 'low')) %>%
  group_by(sample, method, mapping.status) %>%
  summarize(n=sum(n)) %>% 
  group_by(sample, method) %>% mutate(prop=n/sum(n)) %>%
  filter(mapping.status=='high')
gp.df %>% group_by(method) %>% summarize(prop=mean(prop))
ggp.mapq = ggplot(gp.df, aes(x=sample, y=prop, fill=method)) +
  geom_bar(stat='identity', position=position_dodge(.8), width=.7) + 
  scale_fill_manual(values=col.pal, name='genome-mapper') + 
  theme_bw() + coord_flip() + 
  ylab('proportion of reads with mapping quality greater than 0') + 
  ylim(0,1)
ggp.mapq

## sanity check: total number of reads is the same
map.df %>% group_by(sample, method) %>% summarize(n=sum(n))

library(grid)
library(gtable)
library(egg)

pdf(paste0('human-mapping-stats', out.pdf, '.pdf'), 8, 5)
titcap.l =  theme(plot.title=element_text(hjust=-.04, vjust=-3), axis.title.y=element_blank())
grid.arrange(ggp.perfect + ggtitle('A)') + titcap.l + guides(fill='none'),
             ggp.mapq + ggtitle('B)') + titcap.l + theme(legend.position='bottom') +
             guides(fill=guide_legend(nrow=2)),
             heights=c(5,7))
dev.off()

map.df %>% group_by(sample, method) %>%
  summarize(unmapped=n[mapq==-1]/sum(n),
            align.perfect=sum(n[perfect])/sum(n),
            mapq.gt.0=sum(n[mapq>0])/sum(n),
            .groups='drop') %>%
  select(-sample) %>% group_by(method) %>% summarize_all(mean) %>%
  kable(digits=4)
