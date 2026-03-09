#!/bin/bash
#SBATCH -N 1
#SBATCH -n 32
#SBATCH -t 72:00:00
#SBATCH --mem=64GB

# CONFIGURATION:

INPUT_DIR="/route/to/aligned/sam/files/directory"
SAMPLES_LIST_DIR="/route/to/samples/list/directory"
OUTPUT_DIR="/output/directory"

# MODULE LOADING:

module load samtools/1.19

# EXECUTION:

while read -r sample; do
    [[ -z "$sample" || "$sample" =~ ^# ]] && continue
    samtools view -b "${INPUT_DIR}/${sample}.sam" -o "${OUTPUT_DIR}/${sample}.bam" --threads 31
    samtools sort -@ 31 "${OUTPUT_DIR}/${sample}.bam" -o "${OUTPUT_DIR}/${sample}.bam"  # overwritting
    samtools stats -@ 31 "${OUTPUT_DIR}/${sample}.bam" | grep ^SN | cut -f 2- > "${OUTPUT_DIR}/stats/${sample}.txt"  # saving statistics in 'stats' folder
done < "${SAMPLES_LIST_DIR}/samples_list.txt"