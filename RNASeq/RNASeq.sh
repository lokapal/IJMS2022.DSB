#!/bin/sh
# script to get and process H.sapiens HEK293T hg38 gene expression in TPM
# (C) Yuri Kravatsky, lokapal@gmail.com, jiri@eimb.ru
# Input:  H.sapiens HEK293T RNA-Seq data: GEO GSE130262: SRR8950022, SRR8950023
# Output: 1. HEK293.RNAseq1.genes.results    HEK293T hg38 gene expression replica 1
#            HEK293.RNAseq2.genes.results    HEK293T hg38 gene expression replica 1
#         2. HEK293.hg38.TPM                 joined replicas expression all and average values
#            HEK293.hg38.mean.TPM            HEK293T hg38 gene expression only average values per gene
#
# Dependency tools:
# 1. NCBI SRA Toolkit     https://github.com/ncbi/sra-tools
# 2. Trimmomatic 0.36     http://www.usadellab.org/cms/?page=trimmomatic
# 3. RSEM 1.3.1           https://github.com/deweylab/RSEM
# 4. STAR 2.6.1c          https://github.com/alexdobin/STAR
# 5. R with libraries     tibble, dplyr
# Get data from NCBI SRA Archive
fastq-dump --split-files --gzip SRR8950022
fastq-dump --split-files --gzip SRR8950023
mv SRR8950022_1.fastq.gz HEK293.rep1.fastq.gz
mv SRR8950023_1.fastq.gz HEK293.rep2.fastq.gz
# Quality trimming
trimmomatic SE -threads 20 *.rep1.fastq.gz rep1.filtered.fastq.gz LEADING:18 TRAILING:18 SLIDINGWINDOW:4:22 MINLEN:20
trimmomatic SE -threads 20 *.rep2.fastq.gz rep2.filtered.fastq.gz LEADING:18 TRAILING:18 SLIDINGWINDOW:4:22 MINLEN:20
# Expression calculation per replicate
rsem-calculate-expression --star -p 20 --output-genome-bam --calc-ci --ci-memory 30720 --star-gzipped-read-file rep1.filtered.fastq.gz /usr/local/genomes/RNASeq/ref/hg38 HEK293.RNAseq1
rm -f *.bam
rsem-calculate-expression --star -p 20 --output-genome-bam --calc-ci --ci-memory 30720 --star-gzipped-read-file rep2.filtered.fastq.gz /usr/local/genomes/RNASeq/ref/hg38 HEK293.RNAseq2
rm -f *.bam
# Calculate average expression values between replicas
./lib/mean_expr.R
