#!/bin/sh
# script to find intersections between repicates, to filter DFAM, and to assign intersections to genes
# and their expression values 
#
# (C) Yuri Kravatsky, lokapal@gmail.com, jiri@eimb.ru
#
# Input:  1. table_hg38.rep1.txt replicate 1 table (output of the script 02.mapping.sh)
#         2. table_hg38.rep2.txt replicate 2 table (output of the script 02.mapping.sh changed to process rep2)
#         3. lib/HEK293.hg38.mean.TPM gene expression values from RNAseq output (e.g. output from RNASeq directory) 
# Output: 1. HEK293.DSB.hg38.nodfam.intersect.txt        HEK293T DSB genome mappings intersections rep1 and rep2, DFAM removed table
#         2. HEK293.DSB.hg38.nodfam.intersect.bedGraph   HEK293T DSB hg38 profile for genome browsers
#         3. DSB.HEK293.hg38.intersect.all.txt           HEK293T DSB intersect noDFAM with genes and expression added, all entries (with genes intersections and empty)
#         4. DSB.HEK293.hg38.intersect.genes.txt         HEK293T DSB intersect noDFAM with genes and expression added, genes intersections only
#
# Dependency tools:
# 1. bedtools 2.29.1      https://bedtools.readthedocs.io/en/latest/
# 2. BEDOPS 2.4.40        https://bedops.readthedocs.io/en/latest/
# 3. R with libraries refGenome and dplyr
# 4. DFAM database hg38   https://www.dfam.org/releases/Dfam_3.4/annotations/hg38/hg38_dfam.nrph.hits.gz
# 5. hg38 genome annotation in GTF format, unpacked, classic chromosome names like "chr1, chr2, ..."
#                         http://ftp.ensembl.org/pub/current_gtf/homo_sapiens/
#
# Part 1. Find intersection between replicates, remove mappings that are completely inside DFAM entries
# get DFAM database and convert it to bed format
wget https://www.dfam.org/releases/current/annotations/hg38/hg38.nrph.hits.gz
./lib/dfam2bed.pl hg38.nrph.hits.gz
mv hg38.nrph.hits.bed  hg38_dfam.bed
sort -k1,1 -k2,2n -k3,3n hg38_dfam.bed -o hg38_dfam.bed
#
mv table_hg38.rep1.txt rep1.txt
mv table_hg38.rep2.txt rep2.txt
./lib/table2bedg.pl rep1.txt
./lib/table2bedg.pl rep2.txt
sort -k1,1 -k2,2n -k3,3n rep1.bedGraph -o rep1.bedGraph
sort -k1,1 -k2,2n -k3,3n rep2.bedGraph -o rep2.bedGraph
bedtools intersect -sorted -v -f 1.0 -a rep1.bedGraph -b hg38_dfam.bed  > DSB_hg38_nodfam.rep1.bedGraph
bedtools intersect -sorted -v -f 1.0 -a rep2.bedGraph -b hg38_dfam.bed  > DSB_hg38_nodfam.rep2.bedGraph
bedtools intersect -a DSB_hg38_nodfam.rep1.bedGraph -b DSB_hg38_nodfam.rep2.bedGraph -wb > intersect_reps.txt
./lib/overlist2bed_mean.pl intersect_reps.txt
sort -k1,1 -k2,2n -k3,3n intersect_reps.bedGraph -o intersect_reps.bedGraph
sort -k1,1 -k2,2n -k3,3n intersect_reps.bed -o intersect_reps.bed
bedtools merge -c 4 -o mean -i intersect_reps.bed > intersect_merged.bed
./lib/mergebedgraphmean.sh intersect_reps.bedGraph HEK293.DSB.hg38.nodfam.intersect.bedGraph
./lib/addSubseq.pl /usr/local/genomes/hg38.mfa intersect_merged.bed > HEK293.DSB.hg38.nodfam.intersect.txt
sort -k1,1 -k2,2n -k3,3n HEK293.DSB.hg38.nodfam.intersect.txt -o HEK293.DSB.hg38.nodfam.intersect.txt
#
# Part 2. Assign genes to HEK293T DSB-associated intersected replicates noDFAM mappings
cp HEK293.DSB.hg38.nodfam.intersect.txt DSB.bed
# Annotation h38.gtf should be in /usr/local/genomes
Rscript ./lib/getgenes.R
sort -k1,1 -k2,2n -k3,3n hg38.genes.bed -o hg38.genes.bed
bedtools intersect -wa -wb -a DSB.bed -b hg38.genes.bed > DSB.genetable
bedtools intersect -v -wa -a DSB.bed -b hg38.genes.bed > DSB.empty
./lib/tbl_addexpr.pl DSB.genetable ./lib/HEK293.hg38.mean.TPM
cat DSB.genetable.expr DSB.empty > DSB.HEK293.hg38.intersect.all.txt
mv DSB.genetable.expr DSB.HEK293.hg38.intersect.genes.txt
sort -k1,1V -k2,2n -k3,3n DSB.HEK293.hg38.intersect.genes.txt -o DSB.HEK293.hg38.intersect.genes.txt
sort -k1,1V -k2,2n -k3,3n DSB.HEK293.hg38.intersect.all.txt -o DSB.HEK293.hg38.intersect.all.txt
