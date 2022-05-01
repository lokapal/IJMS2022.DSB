#!/bin/sh
# script to get raw reads and remove 4C-rDNA primers/adapters from SE NGS reads
# (C) Yuri Kravatsky, lokapal@gmail.com, jiri@eimb.ru
# Input:  GEO GSE121413: GSM3434713, GSM3434714
# Output: 1. final.rep1.fastq.gz   rep 1 GSM3434713 reads ready for alignment
#            final.rep2.fastq.gz   rep 2 GSM3434714 reads ready for alignment
#
# Dependency tools:
# 1. NCBI SRA Toolkit     https://github.com/ncbi/sra-tools
# 2. cutadapt  3.5        https://cutadapt.readthedocs.io/en/stable/
#
# Get data from SRA Archive
# SRX4900048:GSM3434713 4C-rDNA HEK293T replicate 1
fastq-dump --split-files --gzip SRR8072070
# SRX4900049:GSM3434714 4C-rDNA HEK293T replicate 2
fastq-dump --split-files --gzip SRR8072071
mv SRR8072070_1.fastq.gz rep1.fastq.gz
mv SRR8072071_1.fastq.gz rep2.fastq.gz
#
# cut adapters/primers from replicate 1
cutadapt -j 20 --trim-n --times=4 --minimum-length 20 -q 26 --untrimmed-output=4C_notail.fastq.gz \
-g GCCTAAGCCTGCTGAGAACTTTC -g CAGCATTCTGTAGGGAGATCAAATC -a GAAAGTTCTCAGCAGGCTTAGGC -a GATTTGATCTCCCTACAGAATGCTG \
-o 4C_tail.fastq.gz rep1.fastq.gz
cutadapt -j 20 --trim-n --times=4 --minimum-length 20 -q 26 --untrimmed-output=cutdebug.fastq.gz \
-g TCTTTGAAAAAAATCCCAGAAGTGGT -g AAGTCCAGAAATCAACTCGCCAGT -a ACTGGCGAGTTGATTTCTGGACTT -a ACCACTTCTGGGATTTTTTTCAAAGA \
-o 4C_head.fastq.gz 4C_notail.fastq.gz
cat 4C_tail.fastq.gz 4C_head.fastq.gz > final.rep1.fastq.gz
# cut adapters/primers from replicate 2
cutadapt -j 20 --trim-n --times=4 --minimum-length 20 -q 26 --untrimmed-output=4C_notail.fastq.gz \
-g GCCTAAGCCTGCTGAGAACTTTC -g CAGCATTCTGTAGGGAGATCAAATC -a GAAAGTTCTCAGCAGGCTTAGGC -a GATTTGATCTCCCTACAGAATGCTG \
-o 4C_tail.fastq.gz rep2.fastq.gz
cutadapt -j 20 --trim-n --times=4 --minimum-length 20 -q 26 --untrimmed-output=cutdebug.fastq.gz \
-g TCTTTGAAAAAAATCCCAGAAGTGGT -g AAGTCCAGAAATCAACTCGCCAGT -a ACTGGCGAGTTGATTTCTGGACTT -a ACCACTTCTGGGATTTTTTTCAAAGA \
-o 4C_head.fastq.gz 4C_notail.fastq.gz
cat 4C_tail.fastq.gz 4C_head.fastq.gz > final.rep2.fastq.gz
# remove leftovers
rm -f cutdebug* 4C*.fastq.gz
