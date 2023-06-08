library(dplyr)
library(ggplot2)
library(knitr)
library(egg)

winsor <- function(x, u){
  if(any(x>u)) x[x>u] = u
  x
}

## dna-brnn annotation for each assembly
brnn = read.table('hprc-dnabrnn.bed', as.is=TRUE, sep='\t', comment.char='')
## sum the amount of bases annotated for each contig
brnn = brnn %>% mutate(contig=gsub('(.*)#(.*)#(.*)', "\\1.\\2\\|\\3", V1)) %>%
  group_by(contig) %>% summarize(centsat=sum(V3-V2+1))

## segdup annotation for each assembly
segdup = read.table('hprc-segdup.bed', as.is=TRUE, sep='\t', comment.char='')
## sum the amount of bases annotated for each contig
segdup = segdup %>% mutate(contig=gsub('(.*)#(.*)#(.*)', "\\1.\\2\\|\\3", V1)) %>%
  group_by(contig) %>% summarize(segdup=sum(V3-V2+1))

## function to parse the log files
readLog <- function(filen){
  log.f = gzfile(filen, "r")
  log = scan(log.f, "", sep="\n")
  close(log.f)
  ## skip uninformative lines
  log = grep("Reference contig mappings", log, value=TRUE, invert=TRUE)
  ## grab contig name or force to NA otherwise
  contigs = gsub('.*id=([^ ]+) .*', '\\1', log)
  contigs = ifelse(contigs == log, NA, contigs)
  ## loop over lines and assign contig to NAs
  cur.contigs = NA
  for(ii in 1:length(contigs)){
    if(is.na(contigs[ii])){
      contigs[ii] = cur.contig
    } else {
      cur.contig=contigs[ii]
    }
  }
  ## get the ambiguous chromosome assignment for contigs that have a unique factor below "infinity"
  uf.map = grep("^  ", log)
  uf.map = tibble(contig=contigs[uf.map],
                  chrom=gsub(' *(.*): .*', '\\1', log[uf.map]),
                  cov=as.numeric(gsub('.*: (.*)', '\\1', log[uf.map])))
  ## get the main contig info: assigned or not, length, coverage, uniqueness factor
  df = grep("^  ", log, invert=TRUE)
  df = tibble(contig=contigs[df],
              status=ifelse(grepl('Assigned', log[df]), 'assigned', 'filtered'),
              assignment=gsub('Assigned contig to (.*): id.*', '\\1', log[df]),
              length=as.numeric(gsub('.*len=([^ ]+) .*', '\\1', log[df])),
              cov=as.numeric(gsub('.*cov=([^ ]+) .*', '\\1', log[df])),
              cov.th=as.numeric(gsub('.*cov=[^ ]+ \\(vs ([^ ]+)\\).*', '\\1', log[df])),
              uf=as.numeric(gsub('.*uf= ?([^ ]+) .*', '\\1', log[df])),
              uf.th=as.numeric(gsub('.*uf= ?[^ ]+ \\(vs ([^ ]+)\\).*', '\\1', log[df])),
              amb=gsub('(*)\\|.+', '\\1', contig)) %>%
    mutate(assignment=ifelse(status=='filtered', NA, assignment))
  ## compute the amount of dna-brnn-ed sequence per ambiguous contig
  df = df %>% merge(brnn, all.x=TRUE) %>%
    mutate(centsat=ifelse(is.na(centsat), 0, centsat),
           prop.centsat=centsat/length)
  ## compute the amount of segdup sequence per ambiguous contig
  df = df %>% merge(segdup, all.x=TRUE) %>%
    mutate(segdup=ifelse(is.na(segdup), 0, segdup),
           prop.segdup=segdup/length)
  ## assign the top chromosome for those with low UF
  uf.ass = uf.map %>% group_by(contig) %>% arrange(desc(cov)) %>% summarize(assignment=head(chrom, 1))
  uf.ass.v = uf.ass$assignment
  names(uf.ass.v) = uf.ass$contig
  chrs.order = c(paste0('chr', c(1:22, 'X', 'Y', 'M', 'Other')), 'none')
  df = df %>% mutate(assignment=ifelse(is.na(assignment), uf.ass.v[as.character(contig)],
                                       assignment),
                     assignment=ifelse(is.na(assignment), 'none', assignment),
                     assignment=ifelse(!(assignment %in% chrs.order), 'chrOther', assignment),
                     assignment=factor(assignment, chrs.order),
                     reason=ifelse(cov<cov.th, 'low-cov', NA),
                     reason=ifelse(uf<uf.th, 'low-uf', reason),
                     reason=ifelse(uf<uf.th & cov<cov.th, 'low-cov-uf', reason),
                     reason=ifelse(status=='assigned', 'assigned', reason),
                     reason=factor(reason,
                                   levels=c('assigned', 'low-cov-uf', 'low-cov', 'low-uf'),
                                   labels=c('assigned', 'low coverage and ambiguous', 'low coverage', 'ambiguous')))
  return(list(df=df, uf.map=uf.map))
}

## read GRCh38 and CHM13 logs
gr.l = readLog("GRCh38-f1g-90-mc-jun1.minigraph.split.log.gz")
ch.l = readLog("chm13-f1g-90-mc-jun1.minigraph.split.log.gz")

## combine them
df = rbind(gr.l$df %>% mutate(graph='GRCh38-based pangenome'),
           ch.l$df %>% mutate(graph='CHM13-based pangenome')) %>%
  mutate(graph=factor(graph, unique(graph)))
uf.map = rbind(gr.l$uf.map %>% mutate(graph='GRCh38-based pangenome'),
               ch.l$uf.map %>% mutate(graph='CHM13-based pangenome')) %>%
  mutate(graph=factor(graph, unique(graph)))

## how many ambigous contigs?
df %>% group_by(graph, status) %>% summarize(n=n(), mean.length=mean(length), total.Gbp=sum(length/1e6))

## why were contigs filtered?
df %>% filter(status=='filtered') %>%
  group_by(graph, reason) %>% summarize(n=n())

## count contigs by length class
df.l = df %>% 
  mutate(length.c=cut(length, c(0,1e5,1e6,Inf), c("<100Kbp", '100Kbp-1Mbp', '>1Mbp'))) %>%
  group_by(graph, reason, length.c) %>% summarize(n=n(), mb=sum(length/1e6))

## bar plots for all contigs
ggp.lenn = ggplot(df.l, aes(x=length.c, y=n, fill=reason)) +
  geom_col() + theme_bw() +
  facet_grid(.~graph) + 
  xlab('contig length') + ylab('number of contigs') +
  scale_fill_brewer(palette='Dark2') + 
  theme(legend.title=element_blank())
ggp.lenn

## bar plot for just the filtered contigs
ggp.lenn.f = df.l %>% filter(reason!='assigned') %>%
  ggplot(aes(x=length.c, y=n, fill=reason)) +
  geom_col() + theme_bw() +
  xlab('contig length') + ylab('number of filtered contigs') +
  facet_grid(.~graph) + 
  scale_fill_brewer(palette='Dark2') + 
  theme(legend.title=element_blank())
ggp.lenn.f

## amount of sequence
ggp.lengb = ggplot(df.l, aes(x=length.c, y=mb/1000, fill=reason)) +
  geom_col() + theme_bw() +
  xlab('contig length') + ylab('number of bases (G)') +
  facet_grid(.~graph) + 
  scale_fill_brewer(palette='Dark2') + 
  theme(legend.title=element_blank())
ggp.lengb

## same but log-scaled
ggp.lengb.log = ggplot(df.l, aes(x=length.c, y=mb*1e6, fill=reason)) +
  geom_col(position=position_dodge(preserve='single')) + theme_bw() +
  xlab('contig length') + ylab('number of bases') +
  scale_y_log10(breaks=10^(seq(1,20,2))) + 
  scale_fill_brewer(palette='Dark2') + 
  facet_grid(.~graph) + 
  theme(legend.title=element_blank())
ggp.lengb.log

## how much of each contig is annotated as centromeric satellite
ggp.centsat = ggplot(df, aes(x=prop.centsat)) +
  geom_histogram() + theme_bw() +
  facet_grid(reason~graph, scales='free') +
  theme(strip.text.y=element_text(angle=0)) +
  xlab('proportion of centromeric satellites') +
  ylab('number of contigs') +
  scale_x_continuous(breaks=seq(0,1,.2))
ggp.centsat

## how much of each contig is annotated as segmental duplication
ggp.segdup = ggplot(df, aes(x=prop.segdup)) +
  geom_histogram() + theme_bw() +
  facet_grid(reason~graph, scales='free') +
  theme(strip.text.y=element_text(angle=0)) +
  xlab('proportion of centromeric satellites') +
  ylab('number of contigs') +
  scale_x_continuous(breaks=seq(0,1,.2))
ggp.segdup

## annotate each contig as "centromeric satellite" or "segmental duplication" (or both)
## if more than 10% of its sequence is annotated as such
ggp.ambchr = df %>% filter(grepl('ambiguous', reason)) %>%
  mutate(graph=factor(graph, levels=levels(graph),
                      labels=gsub(' pangenome', '\npangenome', levels(graph))),
         annotation=ifelse(prop.centsat>0.1, 'centromeric satellite', 'none'),
         annotation=ifelse(prop.segdup>0.1, 'segmental duplication', annotation),
         annotation=ifelse(prop.segdup>0.1 & prop.centsat>0.1, 'both', annotation),
         annotation=factor(annotation, c('centromeric satellite',
                                         'segmental duplication', 'both', 'none'))
         ) %>% 
  ggplot(aes(x=assignment, fill=annotation)) +
  geom_bar() + theme_bw() + 
  facet_grid(graph~.) + 
  xlab('best matching reference chromosome') +
  ylab('number of ambiguous contigs') +
  scale_fill_brewer(name='repeat\nannotation',
                    palette='Set1') + 
  theme(axis.text.x=element_text(angle=45, hjust=1),
        strip.text.y=element_text(angle=0),
        legend.position=c(.01,.99), legend.justification=c(0,1)) +
  guides(fill=guide_legend(direction='horizontal', nrow=2))
ggp.ambchr


pdf('contig.inclusion.pdf', 9, 4)
ggp.lenn.f
ggp.ambchr
ggp.centsat
ggp.segdup
dev.off()
