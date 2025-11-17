#!/bin/bash

# $1 is the assembly file used as the reference genome

# (2) Filter the alignments based on mapping quality (MAPQ ≥ 1) and edit distance (NM < 3)
# The 'filter_bam' utility applies these filters to improve alignment accuracy
#   HiC_$1.bam: input BAM file containing Hi-C read alignments
#   1: minimum mapping quality threshold
#   --nm 3: maximum allowed edit distance (alignments with NM >= 3 are discarded)
#   --threads 14: number of threads to use for parallel processing
# The filtered output is then piped to samtools to produce a new BAM file (HiC_$1_filtered.bam)
#   -b: output in BAM format
#   -@: number of threads used by samtools
/users-d1/silvia.garcia/Programs/HapHiC/utils/filter_bam HiC_$1.bam 1 --nm 3 --threads 14 | samtools view - -b -@ 14 -o HiC_$1_filtered.bam

