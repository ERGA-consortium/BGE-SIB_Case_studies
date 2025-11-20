#!/bin/bash

# $1 is the name of the HiFi reads FASTQ file
# $2 is the name of the Hi-C reads FASTQ prefix (without _1.fq.gz or _2.fq.gz)

# Run hifiasm to assemble the HiFi reads using Hi-C data for phasing
# -o specifies the output prefix
# -t16 uses 16 CPU threads
# --h1 and --h2 specify the paired Hi-C read files
# The command output (both stdout and stderr) is redirected to hic.asm.log
hifiasm -o $1 -t16 --h1 $2_1.fq.gz --h2 $2_2.fq.gz $1.fastq.gz &> hic.asm.log

# Extract contig sequences from the haplotype 1 GFA file and convert them to FASTA format
# Each 'S' line in the GFA represents a sequence: print its name ($2) and sequence ($3)
awk '/^S/{print ">"$2; print $3}' $1.hic.hap1.p_ctg.gfa > $1.hic.hap1.fa

# Extract contig sequences from the haplotype 2 GFA file and convert them to FASTA format
awk '/^S/{print ">"$2; print $3}' $1.hic.hap2.p_ctg.gfa > $1.hic.hap2.fa

# Extract contig sequences from the primary GFA assembly and convert them to FASTA format
awk '/^S/{print ">"$2; print $3}' $1.hic.p_ctg.gfa > $1.hic.fa
