#!/bin/bash
#SBATCH -N 1
#SBATCH -n 16
#SBATCH -t 72:00:00
#SBATCH --mem=16GB

# CONFIGURATION:

INPUT_DIR="/route/to/trimmomatic/filtered/paired/files/directory"
SAMPLES_LIST_DIR="/route/to/samples/list/directory"
OUTPUT_DIR="/output/directory"

# MODULE LOADING:

module load cutadapt/3.5

# EXECUTION:

while read -r sample; do
    forward="${INPUT_DIR}/${sample}.R1_001_p.fq.gz"  # change '.R1_001_p' suffix if necessary
    reverse="${INPUT_DIR}/${sample}.R2_001_p.fq.gz"  # change '.R2_001_p' suffix if necessary
    if [[ -f "$forward" && -f "$reverse" ]]; then
        cutadapt -a TACACTCTTTCCCTACACGACGCTCTTCCGATCT -A GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT -o "${OUTPUT_DIR}/${sample}.1.fq.gz" -p "${OUTPUT_DIR}/${sample}.2.fq.gz" "$forward" "$reverse" -j 16
    fi
done < "${SAMPLES_LIST_DIR}/samples_list.txt"
