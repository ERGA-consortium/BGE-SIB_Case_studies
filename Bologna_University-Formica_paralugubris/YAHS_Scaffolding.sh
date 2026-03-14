#!/bin/bash

# Mapping reads with BWA in HiC mode
bwa-mem2 index Fpar_primary.p_ctg.fa
bwa-mem2 mem -5SPC Fpar_scaffolding_of_primary_assembly/Fpar_primary.p_ctg.fa data_HiC/merged_HiC_data/paralugubris_merged_1.fq.gz data_HiC/merged_HiC_data/paralugubris_merged_2.fq.gz -t 16 > Fpar_scaffolding_of_primary_assembly/Fpar_HiC_mapping.sam
#Filtering the alignment
samtools view -b -o Fpar_HiC_mapping.bam Fpar_HiC_mapping.sam
samtools fixmate -m Fpar_HiC_mapping.bam Fpar_HiC_mapping.fixmate.bam
samtools sort -o Fpar_HiC_mapping.namesort.bam Fpar_HiC_mapping.fixmate.bam
samtools markdup -r Fpar_HiC_mapping.fixmate.bam Fpar_HiC_mapping.markdup.bam
samtools view -b -F 0x904 -o Fpar_HiC_mapping.clean.bam Fpar_HiC_mapping.markdup.bam
samtools sort -n -o Fpar_HiC_mapping.namesorted.bam Fpar_HiC_mapping.clean.bam
#Run YAHS
samtools faidx Fpar_primary.p_ctg.fa
yahs -q 0 Fpar_primary.p_ctg.fa Fpar_HiC_mapping.namesorted.bam
#Evaluation of the scaffolding process with PretextView
bwa-mem2 index yahs.out_scaffolds_final.fa
bwa-mem2 mem -t 16 yahs.out_scaffolds_final.fa ../data_HiC/merged_HiC_data/paralugubris_merged_1.fq.gz ../data_HiC/merged_HiC_data/paralugubris_merged_2.fq.gz > yahs_HiC_aligned.bam
samtools sort -@16 -o yahs_HiC_sorted.bam yahs_HiC_aligned.bam
samtools view -@8 -h yahs_HiC_sorted.bam | PretextMap -o map.pretext --sortby length --sortorder descend --mapq 0
