library(dplyr) ## install with: install.packages('dplyr')
library(ggplot2) ## install with: install.packages('ggplot2')

if(file.exists('freebayes-het-rate.tsv')){
  hets.df = read.table('freebayes-het-rate.tsv', sep='\t', as.is=TRUE, header=TRUE)
} else {
  library(VariantAnnotation) ## install with: BiocManager::install('VariantAnnotation')
  ## install BiocManager with: install.packages('BiocManager')
  chrs = c('dm6.chrX', 'dm6.chr2R', 'dm6.chr3L', 'dm6.chr4',
           'dm6.chr3R', 'dm6.chr2L')

  ## get chromosome lengths
  hdr <- scanVcfHeader('results/surject.16-fruitfly-mc-2022-05-26-d2.vcf.gz')
  seqs = seqlengths(hdr)

  ## surject calls
  het.df = lapply(chrs, function(chrn){
    vcf.params = ScanVcfParam(geno="GT", info=NA, which=GRanges(chrn, IRanges(0, seqs[chrn])))
    vcf = readVcf('results/surject.16-fruitfly-mc-2022-05-26-d2.vcf.gz', param=vcf.params)
    hets = colSums(geno(vcf)$GT == '0/1') + colSums(geno(vcf)$GT == '1/0')
    homs = colSums(geno(vcf)$GT == '1/1')
    tibble(seqnames=chrn, sample=names(hets), het=hets, hom=homs, tot=length(vcf))
  })

  het.surject = bind_rows(het.df)

  ## bwa calls
  chrs = c('chrX', 'chrM', 'chr2R', 'chr3L', 'chrY', 'chr4',
           'chr3R', 'chr2L')
  hdr <- scanVcfHeader('results/bwa.dm6.vcf.gz')
  seqs = seqlengths(hdr)

  het.df = lapply(chrs, function(chrn){
    vcf.params = ScanVcfParam(geno="GT", info=NA, which=GRanges(chrn, IRanges(0, seqs[chrn])))
    vcf = readVcf('results/bwa.dm6.vcf.gz', param=vcf.params)
    hets = colSums(geno(vcf)$GT == '0/1') + colSums(geno(vcf)$GT == '1/0')
    homs = colSums(geno(vcf)$GT == '1/1')
    tibble(seqnames=chrn, sample=names(hets), het=hets, hom=homs, tot=length(vcf))
  })

  het.bwa = bind_rows(het.df)

  ## read 
  hets.df = rbind(
    het.bwa %>% mutate(method='bwa'),
    het.surject %>% mutate(method='surject')
  )

  write.table(hets.df, file='freebayes-het-rate.tsv', sep='\t', row.names=FALSE, quote=FALSE)
}

hets.df = hets.df %>%
  mutate(method=ifelse(grepl('bwa', method), 'BWA on dm6',
                       'vg giraffe on\ncactus pangenome'),
         seqnames=gsub('dm6.', '', seqnames)) %>%
  filter(seqnames!='chrM', seqnames!='chrY') %>%
  mutate(seqnames=factor(seqnames, levels=rev(sort(unique(seqnames)))))
head(hets.df)

pdf('freebayes-het-rate.pdf', 8, 5)

ggplot(hets.df, aes(y=het/(het+hom), x=seqnames, fill=method)) +
  geom_boxplot() +
  scale_fill_brewer(palette='Set1', name='mapping approach') + 
  coord_flip() +
  scale_y_continuous(breaks=seq(0,1,.1)) +
  ylab('proportion of heterozygous calls') +
  xlab('chromosome') + 
  theme_bw() +
  theme(legend.position='bottom')

dev.off()

## average numbers
hets.df %>% group_by(sample, method) %>%
  summarize(het=sum(het), hom=sum(hom)) %>%
  group_by(method) %>% 
  summarize(het.prop=mean(het/(het+hom)))
