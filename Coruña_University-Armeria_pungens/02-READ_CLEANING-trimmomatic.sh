#!/bin/bash
#SBATCH -N 1
#SBATCH -n 32
#SBATCH -t 15:00:00
#SBATCH --mem=50GB

# CONFIGURATION:

INPUT_DIR="/route/to/process_radtags/filtered/files/directory"
SAMPLES_LIST_DIR="/route/to/samples/list/directory"
OUTPUT_DIR="/output/directory"
ADAPT_DIR="/route/to/adapters/fasta/file/directory"

# MODULE LOADING:

module load trimmomatic/0.39

# EXECUTION:

# change '_R1_001.1'/'_R2_001.2' input suffixes if necessary
# change '.R1_001_p'/'.R1_001_up'/'.R2_001_p'/'.R2_001_up' output suffixes if necessary (up = unpaired; p = paired)

while read sample; do
	trimmomatic PE -threads 32 -phred33 "${INPUT_DIR}/${sample}_R1_001.1.fq.gz" "${INPUT_DIR}/${sample}_R2_001.2.fq.gz" "${OUTPUT_DIR}/${sample}.R1_001_p.fq.gz" "${OUTPUT_DIR}/${sample}.R1_001_up.fq.gz" "${OUTPUT_DIR}/${sample}.R2_001_p.fq.gz" "${OUTPUT_DIR}/${sample}.R2_001_up.fq.gz" ILLUMINACLIP:"${ADAPT_DIR}/TruSeq3-PE-2.fa":2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:96
done < "${SAMPLES_LIST_DIR}/samples_list.txt"
