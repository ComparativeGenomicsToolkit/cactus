library(methods)
library(dplyr)

args = commandArgs(TRUE)

## read the 4 sets of variants
b.o = read.table(args[1], as.is=TRUE, sep='\t')
s.o = read.table(args[2], as.is=TRUE, sep='\t')
b.b = read.table(args[3], as.is=TRUE, sep='\t')
s.b = read.table(args[4], as.is=TRUE, sep='\t')
colnames(b.o) = colnames(s.o) = colnames(s.b) = colnames(b.b) = c('chr', 'pos', 'ref', 'alt', 'qual', 'gt', 'gq')

## combine in one data.frame
df = rbind(
  b.o %>% mutate(set='bwa_only'),
  s.o %>% mutate(set='surject_only'),
  b.b %>% mutate(set='bwa_both'),
  s.b %>% mutate(set='surject_both')
)

## set up breaks for the QUAL ranges
bks = c(0, 10^c(-16:10))
bks.ii = 1:(length(bks)-1)
names(bks.ii) = paste0('bk', 1:length(bks.ii))

## library(ggplot2)
## ggplot(df, aes(x=qual, fill=gt=='0/0')) +
##   geom_histogram() +
##   theme_bw() +
##   scale_x_log10() +
##   geom_vline(xintercept=bks)

## format qual scores, split in ranges and count number of variant by gt and set
df.qual = df %>% mutate(call=ifelse(grepl('0', gt), 'het', 'hom alt'),
              call=ifelse(gt=='0/0', 'hom ref', call),
              qual.c=cut(qual, breaks=bks, labels=names(bks.ii),
                         include.lowest=TRUE)) %>%
  group_by(set, call, qual.c) %>% summarize(n=n(), .groups='drop') %>%
  mutate(min.qual=bks[bks.ii[as.character(qual.c)]],
         max.qual=bks[1+bks.ii[as.character(qual.c)]],
         qual.c=paste0(ifelse(min.qual==0, '[', '('), min.qual, ',', max.qual, ']')) %>% 
  ungroup %>% arrange(min.qual) %>%
  mutate(qual.c=factor(qual.c, levels=unique(qual.c)))

bks = seq(0, max(df$gq)+5, 5)
bks.ii = 1:(length(bks)-1)
names(bks.ii) = paste0('bk', 1:length(bks.ii))

## ggplot(df, aes(x=gq, fill=gt=='0/0')) +
##   geom_histogram() +
##   theme_bw() + 
##   geom_vline(xintercept=bks)

## format gq scores, split in ranges and count number of variant by gt and set
df.gq = df %>% mutate(call=ifelse(grepl('0', gt), 'het', 'hom alt'),
                      call=ifelse(gt=='0/0', 'hom ref', call),
                      qual.c=cut(gq, breaks=bks, labels=names(bks.ii),
                         include.lowest=TRUE)) %>%
  group_by(set, call, qual.c) %>% summarize(n=n(), .groups='drop') %>%
  mutate(min.qual=bks[bks.ii[as.character(qual.c)]],
         max.qual=bks[1+bks.ii[as.character(qual.c)]],
         qual.c=paste0(ifelse(min.qual==0, '[', '('), min.qual, ',', max.qual, ']')) %>% 
  ungroup %>% arrange(min.qual) %>%
  mutate(qual.c=factor(qual.c, levels=unique(qual.c)))

df.m = rbind(
  df.qual %>% mutate(qual.type='QUAL'),
  df.gq %>% mutate(qual.type='GQ')
)

## write stats in TSV
write.table(df.m, file=args[5], sep='\t', row.names=FALSE, quote=FALSE)

## find regions specific to each method
library(GenomicRanges)

chr.l = df %>% group_by(chr) %>% summarize(l=max(pos))
chr.l.v = chr.l$l
names(chr.l.v) = chr.l$chr

bin.gr = unlist(tileGenome(chr.l.v, tilewidth=1e4), recursive = TRUE)

gr = GRanges(df$chr, IRanges(df$pos, width=ifelse(is.na(df$ref), 1, nchar(df$ref))))
gr$set = df$set
gr = gr[which(df$gt!='0/0')]

bin.gr$bwa.o = countOverlaps(bin.gr, subset(gr, set=='bwa_only'))
bin.gr$bwa.b = countOverlaps(bin.gr, subset(gr, set=='bwa_both'))
bin.gr$surj.o = countOverlaps(bin.gr, subset(gr, set=='surject_only'))
bin.gr$surj.b = countOverlaps(bin.gr, subset(gr, set=='surject_both'))

spec.gr = c(
  subset(bin.gr, bwa.o == 0 & bwa.b ==0 & surj.o>5),
  subset(bin.gr, surj.o == 0 & surj.b ==0 & bwa.o>5)
)

## write stats in TSV
spec.gr %>% as.data.frame %>% write.table(file=args[6], sep='\t', row.names=FALSE, quote=FALSE)
