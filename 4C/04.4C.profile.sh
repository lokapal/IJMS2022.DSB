#!/bin/sh
# script to create average profile from 2 bam files
# (C) Yuri Kravatsky, lokapal@gmail.com, jiri@eimb.ru
#
# Input:  1. rep1.bam rep2.bam                        alignment files for two replicates
#         2. table_hg38.rep1.txt, table_hg38.rep2.txt table files with results (created as a result of 02.4C.mapping.sh)
# Output: 1. 4C_mean.bedGraph, 4C_HEK293_mean.hg38.bw genome-wide RPKM-normalized hg38 HEK293T 4C-rDNA profile for 
#                                                     genome browsers and/or profile charting
# Dependency tools:
# 1. samtools             http://www.htslib.org
# 2. bedtools 2.29.1      https://bedtools.readthedocs.io/en/latest/
# 3. deepTools            https://deeptools.readthedocs.io/en/develop/
# 4. WiggleTools          https://github.com/Ensembl/WiggleTools
# 5. UCSC Kent's tools    https://hgdownload.cse.ucsc.edu/admin/exe/
#
# find intersections
cp table_hg38.rep1.txt rep1.txt
cp table_hg38.rep2.txt rep2.txt
./lib/table2bedg.pl rep1.txt
./lib/table2bedg.pl rep2.txt
sort -k1,1 -k2,2n -k3,3n rep1.bedGraph -o rep1.bedGraph
sort -k1,1 -k2,2n -k3,3n rep2.bedGraph -o rep2.bedGraph
bedtools intersect -a rep1.bedGraph -b rep2.bedGraph -wa -wb > intersect_reps.txt
# define intersections segments range
./lib/overlist2bed_range.pl intersect_reps.txt
# create BAM files from intersections regions
samtools view -b -L intersect_reps.bed -o rep1_I.bam rep1.bam
samtools view -b -L intersect_reps.bed -o rep2_I.bam rep2.bam
samtools index rep1_I.bam
samtools index rep2_I.bam
# create RPKM-normalized profiles per replicate
bamCoverage --effectiveGenomeSize 2913022398 --skipNAs --normalizeUsing RPKM --ignoreForNormalization chr14 --exactScaling -p 12 -b rep1_I.bam -o rep1_I.bw
bamCoverage --effectiveGenomeSize 2913022398 --skipNAs --normalizeUsing RPKM --ignoreForNormalization chr14 --exactScaling -p 12 -b rep2_I.bam -o rep2_I.bw
# create mean profile, convert to BigWig format
rm -f 4C_mean.bedGraph
wiggletools write_bg 4C_mean.bedGraph mean rep1_I.bw rep2_I.bw
sort -k1,1 -k2,2n -k3,3n 4C_mean.bedGraph -o 4C_mean.bedGraph
bedGraphToBigWig 4C_mean.bedGraph ./lib/hg38.chr.sizes 4C_HEK293_mean.hg38.bw
# remove leftovers
rm -f rep1_I.bam rep2_I.bam intersect_reps.* rep1.txt rep2.txt rep1_I.* rep2_I.* rep1.bedGraph rep2.bedGraph
