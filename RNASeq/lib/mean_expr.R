#!/usr/bin/Rscript

suppressPackageStartupMessages(library(tibble))
suppressPackageStartupMessages(library(dplyr))
namehead = "HEK293.RNAseq"
nametail = ".genes.results"

for (num in 1:2) {
   fname <- paste0(namehead,num,nametail)
   temptbl <- read.table(fname, skip=1, sep = "\t", strip.white = TRUE)
   newname <- paste0("Expr",num)
   colnames(temptbl)[6] <- newname
   if (num == 1) {
         alldata <- as.data.frame (temptbl$newname,row.names=temptbl$V1)
         alldata <- as.data.frame (temptbl[[newname]],row.names=as.vector(temptbl$V1))
         colnames(alldata)[num] <- newname
         next
                 }
   alldata <- add_column(alldata, temptbl[[newname]])
   colnames(alldata)[num] <- newname
                 }

  medvect  <- apply(alldata, 1, median)
#  meanvect <- apply(alldata, 1, mean, trim=0.2) # trimmed mean
  meanvect <- apply(alldata, 1, mean)          # "true" mean

#  collist <- c(2:9)
#  subset <- alldata %>% select(collist)
#write.table(as.data.frame(subset), file="subset.txt", row.names=F, col.names=T, sep="\t", quote=F)

  alldata <- add_column(alldata,medvect,meanvect)

  alldata <- rownames_to_column(alldata, var = "GeneID") %>% as_tibble()

  finalmed  <- alldata %>% select(GeneID,medvect)
  finalmean <- alldata %>% select(GeneID,meanvect)

  write.table(as.data.frame(alldata), file="HEK293.hg38.TPM", row.names=F, col.names=T, sep="\t", quote=F)
  write.table(as.data.frame(finalmed), file="HEK293.hg38.median.TPM", row.names=F, col.names=T, sep="\t", quote=F)
  write.table(as.data.frame(finalmean), file="HEK293.hg38.mean.TPM", row.names=F, col.names=T, sep="\t", quote=F)
