#!/bin/bash
#SBATCH -N 1
#SBATCH -n 32
#SBATCH -t 72:00:00
#SBATCH --mem=64GB

# CONFIGURATION:

INPUT_DIR="/route/to/sorted/bam/files/directory"
POPS_DIR="/route/to/populations/map/directory"
OUTPUT_DIR="/output/directory"

# MODULE LOADING:

module load stacks/2.66

# EXECUTION:

ref_map.pl --samples "$INPUT_DIR" --popmap "${POPS_DIR}/populations.txt" --out-path "$OUTPUT_DIR" -T 32
populations -P "$OUTPUT_DIR" --vcf --popmap "${POPS_DIR}/populations.txt" -t 32  # obtaining SNP file in VCF format