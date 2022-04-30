#!/bin/sh
# script to remove RAFT DSB primers/adapters from paired-end reads
# The main paradigm - to be sure that EcorI/PstI primer/adapter is found in ANY of member of read
# and to keep for further processing only with reads that had this adapter.
# (C) Yuri Kravatsky, lokapal@gmail.com, jiri@eimb.ru
# Input:  GSE201829: GSM6073759, GSM6073760
# Output: 1. final.R1.fastq.gz   R1 reads from PE
#            final.R2.fastq.gz   R2 reads from PE
#
# Dependency tools:
# 1. cutadapt  3.5        https://cutadapt.readthedocs.io/en/stable/
#
#>primer_ecorI_pstI at 5'
#CCGAATTCTCCTTATACTGCAGGGG
#>primer_hindIII_NotI - Sau3A
#CCCAAGCTTAAGCGGCCGCAAAC
# With --pair-filter=both --discard-untrimmed, the pair is discarded if both reads do not contain an adapter
#
# GSE201829: SRX15044220
# GSE201829: SRX15044221
# get data, if required
# fastq-dump --split-files --gzip SRRXXXXXXX
# and rename to meaningful names
# Remove MANDATORY 5' EcorI/PstI adapter from R1 ___OR___ R2. RAFT DSB is adjacent to EcorI/PstI adapter.
cutadapt --cores=22 --trim-n --times=5 --minimum-length=20 --quality-cutoff=26 --discard-untrimmed --pair-filter=both \
-g CCGAATTCTCCTTATACTGCAGGGG -G CCGAATTCTCCTTATACTGCAGGGG \
-o trimmed_ecor.R1.fastq.gz -p trimmed_ecor.R2.fastq.gz DSB_HEK293T.rep1.R1.fastq.gz DSB_HEK293T.rep1.R2.fastq.gz
# Remove all other adapter variants - direct and/or reverse complement
cutadapt --cores=22 --times=12 --minimum-length=20 \
-a CCCAAGCTTAAGCGGCCGCAAACX -a GTTTGCGGCCGCTTAAGCTTGGGX -a CCCCTGCAGTATAAGGAGAATTCGGX \
-g XCCCAAGCTTAAGCGGCCGCAAAC -g XGTTTGCGGCCGCTTAAGCTTGGG -g XCCCCTGCAGTATAAGGAGAATTCGG \
-A CCCAAGCTTAAGCGGCCGCAAACX -A GTTTGCGGCCGCTTAAGCTTGGGX -A CCCCTGCAGTATAAGGAGAATTCGGX \
-G XCCCAAGCTTAAGCGGCCGCAAAC -G XGTTTGCGGCCGCTTAAGCTTGGG -G XCCCCTGCAGTATAAGGAGAATTCGG \
-o trimmed.R1.fastq.gz -p trimmed.R2.fastq.gz trimmed_ecor.R1.fastq.gz trimmed_ecor.R2.fastq.gz
# Remove all incomplete adapters
cutadapt --cores=22 --times=12 --minimum-length=20 \
-g file:lib/5adapters.fa -G file:lib/5adapters.fa \
-a file:lib/3adapters.fa -A file:lib/3adapters.fa \
-o filtered.R1.fastq.gz -p filtered.R2.fastq.gz trimmed.R1.fastq.gz trimmed.R2.fastq.gz
#
cutadapt --cores=22 --times=4 --minimum-length=20 \
-a CCCCTGCAGTATAAGGAGAATTCGG -A CCCCTGCAGTATAAGGAGAATTCGG \
-o final.R1.fastq.gz -p final.R2.fastq.gz filtered.R1.fastq.gz filtered.R2.fastq.gz
# remove leftovers
rm -f trimmed*.fastq.gz filtered*.fastq.gz
