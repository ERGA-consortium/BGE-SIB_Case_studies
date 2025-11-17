#!/bin/bash

# $1 is the assembly file used as the reference genome
# $2 is the prefix of the paired-end Hi-C reads (without _1.fq.gz or _2.fq.gz)

# Index the reference genome for alignment with BWA
bwa index $1

# Align Hi-C reads to the reference genome using BWA-MEM
# -5SP are options recommended for Hi-C data (enables chimeric read detection)
# -t 28 uses 28 CPU threads for faster processing
# samblaster is used to mark or remove duplicate reads
# samtools view converts the SAM output to BAM format, keeping headers (-h)
# -@ 14 uses 14 threads for samtools
# -F 3340 filters out unwanted read types (secondary alignments, unmapped reads, etc.)
# The final output is written to HiC_$1.bam
bwa mem -5SP -t 28 $1 $2_1.fq.gz $2_2.fq.gz | samblaster | samtools view - -@ 14 -S -h -b -F 3340 -o HiC_$1.bam

