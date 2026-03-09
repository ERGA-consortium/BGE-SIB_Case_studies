#!/bin/bash
#SBATCH -N 1
#SBATCH -n 16
#SBATCH -t 72:00:00
#SBATCH --mem=16GB

# CONFIGURATION:

INPUT_DIR="/route/to/cutadapt/filtered/files/directory"
SAMPLES_LIST_DIR="/route/to/samples/list/directory"
OUTPUT_DIR="/output/directory"

# MODULE LOADING:

module fastp/0.22.0

# EXECUTION:

while read -r sample; do
    forward="${INPUT_DIR}/${sample}.1.fq.gz"  # change '.1' suffix if necessary
    reverse="${INPUT_DIR}/${sample}.2.fq.gz"  # change '.2' suffix if necessary
    if [[ -f "$forward" && -f "$reverse" ]]; then
        fastp -i "$forward" -I "$reverse" --trim_poly_x -o "${OUTPUT_DIR}/${sample}.1.fq.gz" -O "${OUTPUT_DIR}/${sample}.2.fq.gz" -w 16
    fi
done < "${SAMPLES_LIST_DIR}/samples_list.txt"
