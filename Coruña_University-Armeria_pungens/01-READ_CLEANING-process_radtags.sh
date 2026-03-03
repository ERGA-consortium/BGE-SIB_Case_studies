#!/bin/bash
#SBATCH -N 1
#SBATCH -n 32
#SBATCH -t 15:00:00
#SBATCH --mem=50GB

# CONFIGURATION:

INPUT_DIR="/route/to/process_radtags/filtered/files/directory"
SAMPLES_LIST_DIR="/route/to/samples/list/directory"
OUTPUT_DIR="/output/directory"

# MODULE LOADING:

module load stacks/2.66

# EXECUTION:

while IFS= read -r sample; do
	forward="${RAW_DATA_DIR}/${sample}_R1_001.fastq.gz"  # change '_R1_001' suffix if necessary
	reverse="${RAW_DATA_DIR}/${sample}_R2_001.fastq.gz"  # change '_R2_001' suffix if necessary
	process_radtags -1 "$forward" -2 "$reverse" -i gzfastq -o "$OUTPUT_DIR" --renz-1 pstI --renz-2 bgIII --threads 32 -r -c -q
done < "${SAMPLES_LIST_DIR}/samples_list.txt"
