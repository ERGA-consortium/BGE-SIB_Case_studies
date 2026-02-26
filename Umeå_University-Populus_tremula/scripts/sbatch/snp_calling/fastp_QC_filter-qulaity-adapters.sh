#!/bin/bash -l

#SBATCH -A hpc2n	
#SBATCH -n 1
#SBATCH -c 2
#SBATCH -t 1-00:00:00
#SBATCH -J fastp
#SBATCH --output=reports/sbatch/fastp/sbatch_R-%x_%j-%a.out
#SBATCH --error=reports/sbatch/fastp/sbatch_R-%x_%j-%a.err
#SBATCH -a 0-39

ml GCC/12.3.0
ml GCC/13.2.0
ml fastp/0.23.4
ml multiqc/1.22.3

data_dir="data/fastq/aspen/NordAsp/raw/reads/"
out_dir="data/fastq/aspen/NordAsp/trimmed/"
report_dir="reports/fastp/"

all_file_ids=($(cat info_files/id_fastq_files.txt))

file_id=$(echo ${all_file_ids[$SLURM_ARRAY_TASK_ID]})
sample_id=$(echo ${file_id} | cut -d"_" -f1)

fastp \
    -i ${data_dir}${file_id}_1.fq.gz \
    -I ${data_dir}${file_id}_2.fq.gz \
    -o ${out_dir}${sample_id}.R1.fq.gz \
    -O ${out_dir}${sample_id}.R2.fq.gz \
    -h ${report_dir}html/${sample_id}.html \
    -j ${report_dir}json/${sample_id}.json

