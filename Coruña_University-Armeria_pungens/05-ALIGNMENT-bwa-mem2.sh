#!/bin/bash
#SBATCH -N 1
#SBATCH -n 32
#SBATCH -t 72:00:00
#SBATCH --mem=100GB

# CONFIGURATION:

INPUT_DIR="/route/to/decompressed/filtered/files/directory"
SAMPLES_LIST_DIR="/route/to/samples/list/directory"
OUTPUT_DIR="/output/directory"
REF_DIR="/route/to/reference/genome/fasta/file/directory"

# MODULE LOADING:

module load bwa-mem2/2.2.1

# EXECUTION:

while read -r sample; do
    forward="${INPUT_DIR}/${sample}.1.fq"  # change '.1' suffix if necessary
    reverse="${INPUT_DIR}/${sample}.2.fq"  # change '.2' suffix if necessary
    out="${OUTPUT_DIR}/${sample}.sam"
    if [[ -f "$forward" && -f "$reverse" ]]; then
        bwa-mem2 mem -t 32 "${REF_DIR}/ena_PRJEB86135_sequence.fasta" "$forward" "$reverse" > "$out"
    fi
done < "${SAMPLES_LIST_DIR}/samples_list.txt"
