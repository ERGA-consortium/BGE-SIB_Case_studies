#!/bin/bash

# --- Variables (Generalize here) ---
THREADS=16
ASSEMBLY="polished_assembly.fasta"
HIC_R1="raw_hic_R1.fastq.gz"
HIC_R2="raw_hic_R2.fastq.gz"
N_CHROMS=5  # Change to your target chromosome number
RESTRICTION_ENZYMES="GATC,GANTC,CTNAG,TTAA"

# 1. Adapter trimming (fastp v0.23.4)
module load fastp/0.23.4-GCC-13.2.0
fastp --detect_adapter_for_pe \
      -i "$HIC_R1" -I "$HIC_R2" \
      -o R1_clean.fastq.gz -O R2_clean.fastq.gz \
      --thread 8

# 2. Mapping and basic filtering (samblaster + samtools)
bwa index $ASSEMBLY

# Creates the initial HiC.bam
bwa mem -5SP -t "$THREADS" "$ASSEMBLY" R1_clean.fastq.gz R2_clean.fastq.gz | \
      samblaster | samtools view -@ 8 -S -h -b -F 3340 -o HiC.bam

# 3. Advanced Filtering (filter_bam utility)
# Generates the specialized HiC.filtered.bam for HapHiC
/home/jmorcill/HapHiC/utils/filter_bam HiC.bam 1 --nm 3 --threads 14 | \
      samtools view -b -@ 14 -o HiC.filtered.bam

# 4. HapHiC Pipeline
haphic pipeline "$ASSEMBLY" HiC.filtered.bam "$N_CHROMS" \
       --RE "$RESTRICTION_ENZYMES" \
      --correct_nrounds 2 --threads "$THREADS"

cd /.../04.build/ 
$ bash juicebox.sh
