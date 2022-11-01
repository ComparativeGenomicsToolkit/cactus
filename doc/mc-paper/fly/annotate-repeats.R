library(sveval) ## install with: BiocManager::install('jmonlong/sveval')
library(GenomicRanges) ## install with: BiocManager::install('GenomicRanges')
library(Biostrings) ## install with: BiocManager::install('Biostrings')
library(dplyr) ## install with: install.packages('dplyr')
## install BiocManager with: install.packages('BiocManager')

species='drosophila'
docker.image='jmonlong/repeatmasker:release-4.0.9-p2'
nb.cores = 16

svs = readRDS('16-fruitfly-mc-2022-05-26-d2.100samples.decomposed.svs.site.rol90.insd100.rds')
svs$id = paste0('sv', 1:length(svs))

svs.df = svs %>% as.data.frame %>%
  filter(size>=40) %>% 
  group_by(svsite, type, clique, seqnames) %>%
  arrange(desc(af)) %>% 
  summarize(start=min(start), end=max(end),
    size.min=min(size), size.max=max(size),
    size=size[1], alleles=n(),
    ac=sum(ac), af=sum(af),
    id=id[1],
    .groups='drop')

svs = subset(svs, id %in% svs.df$id)

seqs = svs$alt
seqs[which(svs$type %in% c("DEL", "INV", "DUP"))] = svs$ref[which(svs$type %in% 
                                                                        c("DEL", "INV", "DUP"))]
names(seqs) = svs$id
seqs = as(seqs, "XStringSet")

temp.fa = paste0(tempfile(), ".fa")
writeXStringSet(seqs, temp.fa, format = "FASTA")

temp.dir = dirname(temp.fa)
system2("docker", c("run", "-t", "-v", paste0(temp.dir, 
                                              ":", temp.dir), docker.image, "RepeatMasker", temp.fa, 
                    "--species", species, "-pa", nb.cores))

rmout = utils::read.table(paste0(temp.fa, ".out"), skip = 3, as.is = TRUE, fill = TRUE, header=FALSE)
rmout = rmout[, c(5:7, 10, 11)]
colnames(rmout) = c("id", "start", "end", "repeat.name", 
                    "repeat.class.family")
head(rmout)

rmout$rm.w = rmout$end - rmout$start
rmout = rmout %>% group_by(id, repeat.name, repeat.class.family) %>%
  summarize(rm.w=sum(rm.w)) %>%
  group_by(id) %>% arrange(desc(rm.w)) %>% 
  do(head(., 1)) %>%
  as.data.frame(stringsAsFactors = FALSE)
rownames(rmout) = rmout$id

svs$rmsk.classfam = rmout[svs$id, "repeat.class.family"]
svs$rmsk.name = rmout[svs$id, "repeat.name"]
svs$rmsk.cov = rmout[svs$id, "rm.w"]/GenomicRanges::width(seqs)

saveRDS(svs, file='16-fruitfly-mc-2022-05-26-d2.100samples.decomposed.svs.site.rol90.insd100.rmsk.rds')
