#!/usr/bin/Rscript
suppressPackageStartupMessages(library(refGenome))
suppressPackageStartupMessages(library(dplyr))
gtf <- ensemblGenome()
basedir(gtf) <- "/usr/local/genomes"
read.gtf(gtf, filename="hg38.gtf")
gpe <- getGenePositions(gtf)
final<-gpe %>% select(seqid, start, end, gene_name, gene_id, strand)
write.table(as.data.frame(final), file="hg38.genes.bed", row.names=F, col.names=F, sep="\t", quote=F)
