#!/bin/sh
# script to align to genome previously filtered from adapters, primers and low quality paired end reads
# (C) Yuri Kravatsky, lokapal@gmail.com, jiri@eimb.ru
# Input:  1. GSE201829 HEK293T DSB RAFT paired end reads after adapter removing
# Output: 1. table_hg38.txt the file with the genome-wide hg38 DSB mappings: chromosome coordinates, alignment length, coverage, reads, sequence
#         2. rep1.bam       sorted alignment file that doesn't contain unaligned reads for first replicate
#            rep1.bam.bai   index file for rep1.bam
#
# Dependency tools:
# 1. bwa                  http://bio-bwa.sourceforge.net/
# 2. samtools             http://www.htslib.org
# 3. tabix from samtools 
# 4. Perl with BIO::DB::Fasta BioPerl library
#
bwa mem -t 20 /usr/local/genomes/hg38.mfa final.R1.fastq.gz final.R2.fastq.gz | samtools view -bS -F 4 -o hits.bam
samtools sort -@ 10 hits.bam -O BAM -o hits_sorted.bam
mv hits_sorted.bam hits.bam
samtools mpileup -f /usr/local/genomes/hg38.mfa hits.bam -o pileup.txt
samtools view -@ 10 -O SAM hits.bam -o hits.sam
mv hits.bam rep1.bam
samtools index -@ 10 rep1.bam
sort -k1,1 -k2,2n pileup.txt | bgzip -c > compressed.pileup.gz
tabix -s 1 -b 2 -e 2 compressed.pileup.gz
perl ./lib/makeTablePileup.pl pileup.txt > table.txt
perl ./lib/addNumOfReads.pl hits.sam table.txt > tableA.txt
perl ./lib/addSubseq.pl /usr/local/genomes/hg38.mfa tableA.txt > table_hg38.rep1.txt
rm -f *.sam compressed.* pileup.*
