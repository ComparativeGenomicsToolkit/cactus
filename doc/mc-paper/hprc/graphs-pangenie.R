library(dplyr)
library(ggplot2)
library(tidyr)
library(RColorBrewer)

##
## Pangenie AF across pangenomes
##

## number of SVs by AF
df = read.csv('pangenie-numbers-af.csv', as.is=TRUE)

pan.lev = c('GRCh38-HGSVC', 'GRCh38-HPRC', 'CHM13-HPRC')
pan.col = brewer.pal(6, 'Set1')[c(2,1,3)]
names(pan.col) = pan.lev

df = df %>% mutate(pangenome=factor(pangenome, pan.lev))

ggp.af = ggplot(df, aes(x=min.AF, y=SV.Total, color=pangenome)) +
  geom_line(alpha=.8) + theme_bw() +
  scale_color_manual(values=pan.col) + 
  scale_y_log10(breaks=10^(1:6), labels=as.integer(10^(1:6))) +
  scale_x_continuous(breaks=seq(0,1,.2)) + 
  xlab("minimum allele frequency") +
  ylab("number of structural variants") +
  theme(legend.title=element_blank(), legend.position=c(.01,.01),
        legend.justification=c(.01,.01))
ggp.af


ggp.af.snv = ggplot(df, aes(x=min.AF, y=SNPs, color=pangenome)) +
  geom_line(alpha=.8) + theme_bw() +
  scale_color_manual(values=pan.col) + 
  scale_y_log10(breaks=10^(1:9), labels=as.integer(10^(1:9))) +
  scale_x_continuous(breaks=seq(0,1,.2)) + 
  xlab("minimum allele frequency") +
  ylab("number of SNPs") +
  theme(legend.title=element_blank(), legend.position=c(.01,.01),
        legend.justification=c(.01,.01))
ggp.af.snv

ggp.af.indel = ggplot(df, aes(x=min.AF, y=Indels, color=pangenome)) +
  geom_line(alpha=.8) + theme_bw() +
  scale_color_manual(values=pan.col) + 
  scale_y_log10(breaks=10^(1:6), labels=as.integer(10^(1:6))) +
  scale_x_continuous(breaks=seq(0,1,.2)) + 
  xlab("minimum allele frequency") +
  ylab("number of indels") +
  theme(legend.title=element_blank(), legend.position=c(.01,.01),
        legend.justification=c(.01,.01))
ggp.af.indel

df.sv = df %>% pivot_longer(cols=c(SV.Deletions, SV.Insertions), names_to='type') %>%
  mutate(type=tolower(gsub('SV.(.*)s', '\\1', type)))

ggp.af.type = ggplot(df.sv, aes(x=min.AF, y=value, linetype=type, color=pangenome)) +
  geom_line(alpha=.7, size=1) + theme_bw() +
  scale_color_manual(values=pan.col) + 
  scale_y_log10(breaks=10^(1:6), labels=as.integer(10^(1:6))) +
  scale_x_continuous(breaks=seq(0,1,.2)) + 
  xlab("minimum allele frequency") +
  ylab("number of structural variants") +
  theme(legend.title=element_blank(), legend.position=c(.01,.01),
        legend.justification=c(.01,.01), legend.box = 'horizontal')
ggp.af.type

## per sample
s.df = lapply(list.files('af-per-samples/'), function(ff){
  dd = read.table(paste0('af-per-samples/', ff), as.is=TRUE)
  dd1 = dd[,c(1,2:7)]
  dd2 = dd[,c(1,(2:7)+6)]
  dd3 = dd[,c(1,(2:7)+12)]
  colnames(dd1) = colnames(dd2) = colnames(dd3) = c('sample', colnames(df)[1:6])
  rbind(dd1 %>% mutate(pangenome='GRCh38-HPRC'),
        dd2 %>% mutate(pangenome='CHM13-HPRC'),
        dd3 %>% mutate(pangenome='GRCh38-HGSVC'))
}) %>% bind_rows

s.df = s.df %>% mutate(SVs=SV.Deletions+SV.Insertions+SV.Other,
                       pangenome=factor(pangenome, pan.lev))

s.df.m = s.df %>% group_by(min.AF, pangenome) %>% summarize(SVs=median(SVs))

ggp.af.s = ggplot(s.df, aes(x=factor(min.AF), y=SVs, fill=pangenome)) +
  geom_violin(scale='width') +
  geom_point(data=s.df.m, position=position_dodge(.9), size=1) +
  geom_point(aes(color=pangenome), data=s.df.m, size=.5, position=position_dodge(.9)) +
  scale_fill_manual(values=pan.col) + 
  scale_color_manual(values=pan.col) + 
  xlab("minimum allele frequency") +
  ylab("number of SVs per sample") +
  theme_bw() +
  theme(legend.title=element_blank(),
        legend.position=c(.99,.99),
        legend.justification=c(1,1))
ggp.af.s

s.t.df = s.df %>% pivot_longer(cols=c(SV.Insertions, SV.Deletions),
                      names_to='type') %>%
  mutate(type=tolower(gsub("SV.(.*)s", '\\1', type)))
s.t.df.m = s.t.df %>% group_by(min.AF, pangenome, type) %>%
  summarize(value=median(value))

svt="deletion"
ggp.af.s.del = ggplot(subset(s.t.df, type==svt),
                      aes(x=factor(min.AF), y=value, fill=pangenome)) +
  geom_violin(scale='width') +
  geom_point(data=subset(s.t.df.m, type==svt),
             position=position_dodge(.9), size=1) +
  geom_point(aes(color=pangenome),
             data=subset(s.t.df.m, type==svt),
             size=.5,
             position=position_dodge(.9)) +
  scale_fill_manual(values=pan.col) + 
  scale_color_manual(values=pan.col) + 
  xlab("minimum allele frequency") +
  ylab(paste0("number of ", tolower(svt), "s per sample")) +
  theme_bw() +
  theme(legend.title=element_blank(),
        legend.position=c(.99,.99),
        legend.justification=c(1,1))
ggp.af.s.del

svt="insertion"
ggp.af.s.ins = ggplot(subset(s.t.df, type==svt),
                      aes(x=factor(min.AF), y=value, fill=pangenome)) +
  geom_violin(scale='width') +
  geom_point(data=subset(s.t.df.m, type==svt),
             position=position_dodge(.9), size=1) +
  geom_point(aes(color=pangenome),
             data=subset(s.t.df.m, type==svt),
             size=.5,
             position=position_dodge(.9)) +
  scale_fill_manual(values=pan.col) + 
  scale_color_manual(values=pan.col) + 
  xlab("minimum allele frequency") +
  ylab(paste("number of ", tolower(svt), "s per sample")) +
  theme_bw() +
  theme(legend.title=element_blank(),
        legend.position=c(.99,.99),
        legend.justification=c(1,1))
ggp.af.s.ins

##
## Leave-one-out experiments
##

pan.lev = c('GRCh38-HGSVC', 'GRCh38-HPRC', 'CHM13-HPRC')
pan.col = brewer.pal(6, 'Set1')[c(2,1,3)]
names(pan.col) = pan.lev

## read the files with the genotype concordances
df = lapply(c('hgsvc', 'hg38', 'chm13'), function(pg){
  lapply(list.files(paste0('concordances-', pg)), function(conc.ff){
    read.table(paste0('concordances-', pg, "/", conc.ff), as.is=TRUE, header=TRUE) %>%
      mutate(allele=gsub("concordance_(.+)_.+.tsv", '\\1', conc.ff),
             type=gsub("concordance_.+_(.+).tsv", '\\1', conc.ff),
             pangenome=pg)    
  }) %>% bind_rows
}) %>% bind_rows

df = df %>% 
  filter(pangenome != 'hgsvc' | type != 'large-complex') %>%
  mutate(pangenome=factor(pangenome, c('hgsvc', 'hg38', 'chm13'), pan.lev),
         type=factor(type, c('snp', 'indels', 'large-deletion',
                             'large-insertion', 'large-complex'),
                     labels=c('SNP', 'indel', 'large deletion',
                             'large insertion', 'complex SV')))

ggp.hprc = df %>% group_by(sample) %>% filter(length(unique(pangenome))>1) %>%
  ungroup %>% arrange(pangenome) %>%
  mutate(pangenome=factor(pangenome, unique(pangenome)),
         sample=reorder(sample, weighted_genotype_concordance, mean)) %>%
  ggplot(aes(x=sample, y=weighted_genotype_concordance,
               color=pangenome)) +
  geom_line(aes(group=paste(pangenome, allele))) +
  theme_bw() +
  scale_color_manual(values=pan.col[2:3]) + 
  facet_grid(allele~type, scales='free') +
  theme(axis.text.x=element_text(angle=45, hjust=1)) +
  ylab('weighted genotype concordance')
ggp.hprc

ggp.hprc.bp = df %>% group_by(sample) %>% filter(length(unique(pangenome))>1) %>%
  ungroup %>% arrange(pangenome) %>%
  mutate(pangenome=factor(pangenome, unique(pangenome))) %>%
  ggplot(aes(x=type, y=weighted_genotype_concordance,
             fill=pangenome)) +
  geom_boxplot() +
  theme_bw() +
  scale_fill_manual(values=pan.col[2:3]) +
  facet_grid(allele~., scales='free') + 
  theme(axis.text.x=element_text(angle=20, hjust=1)) +
  xlab('variant type')
  ylab('weighted genotype concordance')
ggp.hprc.bp

ggp.all.eb = df %>%
  group_by(type, pangenome, allele) %>%
  summarize(wg.med=median(weighted_genotype_concordance),
            wg.uq=max(weighted_genotype_concordance),
            wg.lq=min(weighted_genotype_concordance)) %>%
  ggplot(aes(x=type, y=wg.med,
             ymin=wg.lq, ymax=wg.uq,
             color=pangenome)) +
  geom_linerange(position=position_dodge(.6)) +
  geom_point(position=position_dodge(.6), size=2) +
  theme_bw() +
  scale_color_manual(values=pan.col) +
  facet_wrap(~allele, scales='free') + 
  theme(axis.text.x=element_text(angle=30, hjust=1)) +
  ylab('weighted genotype concordance') +
  xlab('variant type')
ggp.all.eb
