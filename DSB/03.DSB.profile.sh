#!/bin/sh
# script to build intersected mean profile from two BAM files
# (C) Yuri Kravatsky, lokapal@gmail.com, jiri@eimb.ru
# Input:  1. alignment sorted BAM files for HEK293T DSB for earch replicate (rep1.bam, rep2.bam)
#         2. HEK293.DSB.hg38.intersect.bed - the file with intersections in BED format (e.g. output from script 04.DSB.genes.sh)
# Output: 1. HEK293.DSB.fseq.intersect.bw Profile built by F-Seq and only for segments that are intersected.
#
# Dependency tools:
# 1. samtools             http://www.htslib.org
# 2. bedtools 2.29.1      https://bedtools.readthedocs.io/en/latest/
# 3. F-seq 1.84           https://fureylab.web.unc.edu/software/fseq/
# 4. WiggleTools          https://github.com/Ensembl/WiggleTools
# 5. UCSC Kent's tools    https://hgdownload.cse.ucsc.edu/admin/exe/
#
samtools view -@ 10 -b -L HEK293.DSB.hg38.intersect.bed -o rep1_I.bam rep1.bam
samtools view -@ 10 -b -L HEK293.DSB.hg38.intersect.bed -o rep2_I.bam rep2.bam
samtools index -@ 10 rep1_I.bam
samtools index -@ 10 rep2_I.bam
bamToBed -ed -i rep1_I.bam > density.rep1.bed
bamToBed -ed -i rep2_I.bam > density.rep2.bed
fseq density.rep1.bed -o rep1
fseq density.rep2.bed -o rep2
cd rep1
../lib/join.rep1.sh
mv rep1.wig ../
rm -f *.wig
cd ../rep2
../lib/join.rep2.sh
mv rep2.wig ../
rm -f *.wig
cd ..
wiggletools write_bg DSB_mean.bedGraph mean rep1.wig rep2.wig
sort -k1,1 -k2,2n -k3,3n DSB_mean.bedGraph -o DSB_mean.bedGraph
bedGraphToBigWig DSB_mean.bedGraph ~/bin/hg38.chr.sizes HEK293.DSB.fseq.intersect.bw
rm -f DSB_mean.bedGraph density*.bed
