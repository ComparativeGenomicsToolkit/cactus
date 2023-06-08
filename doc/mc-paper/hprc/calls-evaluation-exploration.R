library(dplyr) ## install with: install.packages('dplyr')
library(ggplot2) ## install with: install.packages('ggplot2')
library(tidyr) ## install with: install.packages('tidyr')
library(ggrepel) ## install with: install.packages('ggrepel')
library(grid)
library(gtable)
library(egg)
library(RColorBrewer)
library(colorspace)

## some color palettes
set1.pal = brewer.pal(9, 'Set1')
pair.pal = brewer.pal(8, 'Paired')

## evaluation on GRCh38
df.hg002 = read.table('eval-summary-HG002.tsv', sep='\t', as.is=TRUE, header=TRUE)
df.giab = read.table('eval-summary-GIAB_4_2_1_HG00_1_2_5.tsv', sep='\t', as.is=TRUE, header=TRUE)
df.chm13lift = read.table('eval-summary-GIAB_4_2_1_chm13visible.tsv', sep='\t', as.is=TRUE, header=TRUE)
df.chm13liftcmrg = read.table('eval-summary-CMRG_1_0_chm13visible.tsv', sep='\t', as.is=TRUE, header=TRUE)

df = rbind(
  df.hg002,
  df.chm13lift,
  df.chm13liftcmrg,
  df.giab) %>%
  filter(Filter=='ALL', !duplicated(paste(sample, method, region, truthset, Type)))

## evaluation on CHM13
df.chm13 = read.table('eval-summary-HG002_chm13.tsv', sep='\t', as.is=TRUE, header=TRUE)

## to rename regions later
regions_n = c(wg_noinconsistent='Whole genome',
              nofalsedup_in_chm13='in CHM13 pangenome and no false duplications/collapse')

## to rename methods later
mapper.n = c(
  giraffedv='Giraffe on HPRC-GRCh38',
  giraffe_noclip='Giraffe on full HPRC-GRCh38',
  giraffedv_chm13='Giraffe on HPRC-CHM13',
  giraffedv_chm13_lifted='Giraffe on HPRC-CHM13',
  bwadv='BWAMEM on GRCh38',
  bwadv_chm13_lifted='BWAMEM on CHM13',
  dragen='DragenGRAPH on GRCh38+'
)
caller.n = c(
  giraffedv='DeepVariant',
  giraffe_noclip='DeepVariant',
  giraffedv_chm13='DeepVariant',
  giraffedv_chm13_lifted='DeepVariant-lifted',
  bwadv_chm13_lifted='DeepVariant-lifted',
  bwadv='DeepVariant',
  dragen='Dragen'
)
mapper.n.chm13 = c(
  giraffedv='Giraffe on HPRC-CHM13',
  bwadv='BWAMEM on CHM13'
)

## to rename truthsets
truthset_order = c("GIAB_4_2_1",
                   "GIAB_4_2_1_chm13visible",
                   'GIAB_4_2_1_lifted',
                   '20211005_dipcall_z2k',
                   "CMRG_1_0_chm13visible",
                   "CMRG_1_0",
                   "20211005_dipcall-z2k")
truthset_labels = c("GIAB v4.2.1",
                    "GIAB v4.2.1\nCHM13-visible",
                    'GIAB v4.2.1 confident\nregions lifted',
                    '20211005_dipcall_z2k',
                    "CMRG v1.0\nCHM13-visible",
                    "CMRG v1.0",
                    "20211005 dipcall z2k")

## reformat data frames
df = df %>% filter(!grepl('jointcalled', method)) %>%
  mutate(region=factor(region, levels=names(regions_n), labels=regions_n),
         mapper=factor(mapper.n[method], levels=unique(mapper.n)),
         variant_caller=factor(caller.n[method], levels=unique(caller.n)),
         Total.Error=TRUTH.FN+QUERY.FP) %>%
  filter(!is.na(region)) %>% 
  dplyr::rename(TP.baseline=TRUTH.TOTAL, TP.call=TRUTH.TP, type=Type,
                FN=TRUTH.FN, FP=QUERY.FP, recall=METRIC.Recall,
                precision=METRIC.Precision, F1=METRIC.F1_Score) %>% 
  select(region, sample, truthset, reads, mapper, variant_caller, type, TP.baseline,
         TP.call, FN, FP, recall, precision, F1, Total.Error) %>% 
  arrange(region, sample, truthset, reads, mapper, variant_caller, type)

df.chm13 = df.chm13 %>% filter(Filter=='ALL', !grepl('jointcalled', method)) %>%
  mutate(region=factor(region, levels=names(regions_n), labels=regions_n),
         mapper=factor(mapper.n.chm13[method], levels=unique(mapper.n.chm13)),
         variant_caller=factor(caller.n[method], levels=unique(caller.n)),
         Total.Error=TRUTH.FN+QUERY.FP) %>%
  filter(!is.na(region)) %>% 
  dplyr::rename(TP.baseline=TRUTH.TOTAL, TP.call=TRUTH.TP, type=Type,
                FN=TRUTH.FN, FP=QUERY.FP, recall=METRIC.Recall,
                precision=METRIC.Precision, F1=METRIC.F1_Score) %>% 
  select(region, sample, truthset, reads, mapper, variant_caller, type, TP.baseline,
         TP.call, FN, FP, recall, precision, F1, Total.Error) %>% 
  arrange(region, sample, truthset, reads, mapper, variant_caller, type)

## compute stats for SNPs and indels combined
df.r = df %>%
  mutate(truthset=factor(truthset, levels=truthset_order, labels=truthset_labels)) %>% 
  group_by(sample, truthset, mapper, variant_caller, region) %>% 
  summarize(TP=sum(TP.call), total=sum(TP.baseline),
            FP=sum(FP), FN=sum(FN)) %>%
  mutate(precision=TP/(TP+FP), recall=TP/(TP + FN),
         F1=2*precision * recall/(precision + recall))

df.chm13.r = df.chm13 %>%
  mutate(truthset=factor(truthset, levels=truthset_order, labels=truthset_labels)) %>% 
  group_by(sample, truthset, mapper, variant_caller, region) %>% 
  summarize(TP=sum(TP.call), total=sum(TP.baseline),
            FP=sum(FP), FN=sum(FN)) %>%
  mutate(precision=TP/(TP+FP), recall=TP/(TP + FN),
         F1=2*precision * recall/(precision + recall))

##
## figures
##

## main figure: GRCh38-based pangenome and CHM13-based pangenome projected to GRCh38
comb.pal = c(darken(set1.pal[1], .3), lighten(set1.pal[1], .2),
             set1.pal[2:3])
ggp.hprc.main = rbind(
  df.r %>% filter(!grepl('visible', truthset),
                  !grepl('lifted', variant_caller),
                  !grepl('jun22', mapper),
                  !grepl('full', mapper),
                  !grepl('CHM13', mapper),
                  region=='Whole genome') %>% mutate(exp='Evaluation of the\nGRCh38-based\nHPRC pangenome') %>%
  ungroup %>% select(mapper, truthset, sample, F1, exp),
  df.r %>% filter(!grepl('visible', truthset),
                  !grepl('lifted', variant_caller),
                  !grepl('full', mapper),
                  !grepl('jun22', mapper),
                  region=='in CHM13 pangenome and no false duplications/collapse') %>%
  mutate(exp='Evaluation of the\nCHM13-based\nHPRC pangenome') %>%
  ungroup %>% select(mapper, truthset, sample, F1, exp)) %>%
  mutate(truthset=ifelse(grepl('GIAB', truthset), 'GIAB v4.2.1', 'CMRG v1.0'),
         truthset=factor(truthset, levels=c('GIAB v4.2.1', 'CMRG v1.0')),
         exp=factor(exp, levels=unique(exp))) %>%
  arrange(mapper) %>%
  mutate(mapper=factor(mapper, levels=levels(mapper),
                       labels=paste0(levels(mapper), ' and ',
                                     ifelse(grepl('Dragen', levels(mapper)), 'Dragen', 'DeepVariant')))) %>% 
  ggplot(aes(x=F1, y=sample, color=mapper)) +
  geom_point(aes(group=mapper), position=position_dodge(.4), alpha=.8) +
  scale_colour_manual(values=comb.pal) + 
  facet_grid(exp~truthset, scales='free', space='free_y') +
  theme_bw() + 
  theme(legend.position='bottom',
        legend.title=element_blank(),
        strip.text.y=element_text(angle=0)) +
  guides(color=guide_legend(ncol=3))
print(ggp.hprc.main)

## supp: calling on CHM13 and either lifting back to GRCh38 or eval on CHM13
comb.pal = c(set1.pal[1:2])
ggp.hprc.chm13 = rbind(
  df %>% mutate(truthset=factor(truthset, levels=truthset_order,
                                labels=truthset_labels)) %>%
  filter(grepl('visible', truthset),
         type=='SNP',
         grepl('lifted', variant_caller),
         !grepl('jun22', mapper),
         region=='Whole genome') %>%
  mutate(exp='Calls on CHM13 lifted over\nto GRCh38 (SNP only and\non visible variants)') %>%
  ungroup %>% select(mapper, truthset, sample, F1, exp),
  df.chm13.r %>% filter(region=='Whole genome', !grepl('dipcall', truthset)) %>%
  mutate(exp='CMRG CHM13v1.0') %>%
  ungroup %>% select(mapper, truthset, sample, F1, exp),
  df.chm13.r %>% filter(region=='Whole genome', grepl('dipcall', truthset)) %>%
  mutate(exp='dipcall CHM13v2.0') %>%
  ungroup %>% select(mapper, truthset, sample, F1, exp)) %>%
  mutate(truthset=ifelse(grepl('GIAB', truthset), 'GIAB v4.2.1', 'CMRG v1.0'),
         truthset=factor(truthset, levels=c('GIAB v4.2.1', 'CMRG v1.0')),
         exp=factor(exp, levels=unique(exp))) %>%
  arrange(mapper) %>%
  mutate(mapper=factor(mapper, levels=levels(mapper),
                       labels=paste0(levels(mapper), ' and ',
                                     ifelse(grepl('Dragen', levels(mapper)), 'Dragen', 'DeepVariant')))) %>% 
  ggplot(aes(x=F1, y=sample, color=mapper)) +
  geom_point(aes(group=mapper), position=position_dodge(.4), alpha=.8) +
  scale_colour_manual(values=comb.pal) + 
  facet_grid(exp~truthset, scales='free', space='free_y') +
  theme_bw() + 
  theme(legend.position='bottom',
        legend.title=element_blank(),
        strip.text.y=element_text(angle=0)) +
  guides(color=guide_legend(ncol=3))
print(ggp.hprc.chm13)

## supp: calling using the full pangenome vs frequency-filtered on
comb.pal = c(darken(set1.pal[1], .3), lighten(set1.pal[1], .2),
             set1.pal[2:3])
ggp.hprc.noclip = 
  df.r %>% filter(!grepl('visible', truthset),
                  !grepl('lifted', variant_caller),
                  !grepl('jun22', mapper),
                  !grepl('CHM13', mapper),
                  region=='Whole genome') %>% 
  ungroup %>% select(mapper, truthset, sample, F1) %>%
  mutate(truthset=ifelse(grepl('GIAB', truthset), 'GIAB v4.2.1', 'CMRG v1.0'),
         truthset=factor(truthset, levels=c('GIAB v4.2.1', 'CMRG v1.0'))) %>%
  arrange(mapper) %>%
  mutate(mapper=factor(mapper, levels=levels(mapper),
                       labels=paste0(levels(mapper), ' and ',
                                     ifelse(grepl('Dragen', levels(mapper)), 'Dragen', 'DeepVariant')))) %>% 
  ggplot(aes(x=F1, y=sample, color=mapper)) +
  geom_point(aes(group=mapper), position=position_dodge(.4), alpha=.8) +
  scale_colour_manual(values=comb.pal) + 
  facet_grid(.~truthset, scales='free', space='free_y') +
  theme_bw() + 
  theme(legend.position='bottom',
        legend.title=element_blank(),
        strip.text.y=element_text(angle=0)) +
  guides(color=guide_legend(ncol=3))
print(ggp.hprc.noclip)

## number of errors
comb.pal = c(darken(set1.pal[1], .3), lighten(set1.pal[1], .2),
             set1.pal[2:3])
ggp.error.hprc.main = rbind(
  df.r %>% filter(!grepl('visible', truthset),
                  !grepl('lifted', variant_caller),
                  !grepl('jun22', mapper),
                  !grepl('default', mapper),
                  !grepl('CHM13', mapper),
                  region=='Whole genome') %>% mutate(exp='Evaluation of the\nGRCh38-based\nHPRC pangenome') %>%
  ungroup %>% select(mapper, truthset, sample, FP, FN, TP, F1, exp),
  df.r %>% filter(!grepl('visible', truthset),
                  !grepl('lifted', variant_caller),
                  !grepl('default', mapper),
                  !grepl('jun22', mapper),
                  region=='in CHM13 pangenome and no false duplications/collapse') %>%
  mutate(exp='Evaluation of the\nCHM13-based\nHPRC pangenome') %>%
  ungroup %>% select(mapper, truthset, sample, FP, FN, TP, F1, exp)) %>%
  mutate(truthset=ifelse(grepl('GIAB', truthset), 'GIAB v4.2.1', 'CMRG v1.0'),
         truthset=factor(truthset, levels=c('GIAB v4.2.1', 'CMRG v1.0')),
         exp=factor(exp, levels=unique(exp))) %>%
  arrange(mapper) %>%
  mutate(mapper=factor(mapper, levels=levels(mapper),
                       labels=paste0(levels(mapper), ' and ',
                                     ifelse(grepl('Dragen', levels(mapper)), 'Dragen', 'DeepVariant')))) %>%
  select(sample, mapper, exp, truthset, FN, FP) %>%
  pivot_longer(cols=c(FN,FP)) %>%
  group_by(sample, mapper, exp, truthset) %>%
  mutate(ys=ifelse(name=='FN', 0, value[name=='FN']),
         ye=ifelse(name=='FN', value, value+value[name=='FN'])) %>%
  ggplot(aes(x=sample, ymin=ys, ymax=ye, alpha=name, color=mapper)) +
  geom_linerange(aes(group=mapper), position=position_dodge(.8), size=1.3) +
  scale_color_manual(values=comb.pal) +
  scale_alpha_manual(values=c(.5,1), labels=c('false-negative', 'false-positive')) + 
  facet_grid(exp~truthset, scales='free') +
  theme_bw() + coord_flip() +
  ylab('number of erroneous variant calls') + 
  theme(legend.position='bottom',
        legend.title=element_blank(),
        strip.text.y=element_text(angle=0)) +
  guides(color=guide_legend(ncol=3, override.aes=list(size=5)),
         alpha=guide_legend(ncol=1, override.aes=list(size=5)))
print(ggp.error.hprc.main)

pdf('hprc_smallvariant-giab.pdf', 8, 3)
print(ggp.hprc.main)
print(ggp.hprc.chm13)
print(ggp.hprc.noclip)
print(ggp.error.hprc.main)
dev.off()


##
## ROC
##

df = read.table('eval-roc-summary-HG002.tsv', sep='\t', as.is=TRUE, header=TRUE)

df = df %>% filter(!grepl('jointcalled', method)) %>%
  dplyr::rename(type=Type) %>% 
  mutate(region=factor(region, levels=names(regions_n), labels=regions_n),
         mapper=factor(mapper.n[method], levels=unique(mapper.n)),
         variant_caller=factor(caller.n[method], levels=unique(caller.n)),
         truthset=factor(truthset, levels=truthset_order, labels=truthset_labels)
         ) %>% 
  filter(!is.na(region))

methods_colors = brewer.pal(4, 'Set1')[c(1,1:3)]
methods_colors[1] = darken(methods_colors[1], .3)
methods_colors[2] = lighten(methods_colors[2], .2)

df.r = df %>%
  group_by(sample, truthset, mapper, variant_caller, QQ, region) %>% 
  summarize(TP=sum(TRUTH.TP), total=sum(TRUTH.TOTAL),
            FP=sum(QUERY.FP), FN=sum(TRUTH.FN)) %>%
  mutate(precision=TP/(TP+FP), recall=TP/(TP + FN),
         f1=2*precision * recall/(precision + recall))
df.s = df.r %>% 
  filter(grepl('in CHM13', region),
         !grepl('jun22', mapper),
         !grepl('full', mapper),
         sample=='HG002',
         grepl('CMRG', truthset)) %>%
  arrange(mapper, variant_caller) %>% 
  mutate(method=paste(mapper, variant_caller, sep='\n'),
         method=factor(method, unique(method))) %>%
  group_by(method) %>% arrange(desc(f1)) %>%
  mutate(label=ifelse(1:n()==1, paste0(method, '\nF1: ', round(f1, 4)), ''))
lim.expand = .01
pr.min = df.s %>% group_by(method) %>%
  summarize(pr.min=min(max(precision, na.rm=TRUE),
                       max(recall, na.rm=TRUE)), .groups='drop') %>%
  .$pr.min %>% min
df.gpp = df.s %>%
  filter(precision>pr.min-.1*(1-pr.min),
         recall>pr.min-.2*(1-pr.min)) %>% 
  arrange(method, QQ)
xlims = c(min(df.gpp$precision), max(df.gpp$precision))
xlims[2] = xlims[2] + lim.expand*diff(xlims)
ylims = c(min(df.gpp$recall), max(df.gpp$recall))
ylims[2] = ylims[2] + lim.expand*diff(ylims)

pdf('hprc_smallvariant-giab-roc.pdf', 5, 5)
ggp.hprc.roc = ggplot(df.gpp, aes(x=precision, y=recall, color=method)) +
  geom_path(size=1, alpha=.8) +
  theme_bw() +
  theme(legend.position='bottom') +
  scale_color_manual(values=methods_colors, name='method') +
  scale_linetype_manual(name='mapping', values=c(1,2,3)) +
  scale_y_continuous(expand=expansion(mult=.1)) + 
  geom_label_repel(aes(label=label),
                   box.padding = 0.2,
                   label.padding = 0.2,
                   nudge_x=-.005,
                   max.overlaps=10,
                   force=1,
                   force_pull=1,
                   alpha=1,
                   segment.linetype=2,
                   size=3, show.legend=FALSE) +
  xlim(xlims[1], xlims[2]) + 
  ylim(ylims[1], ylims[2]) + 
  guides(color=FALSE)
ggp.hprc.roc
dev.off()

