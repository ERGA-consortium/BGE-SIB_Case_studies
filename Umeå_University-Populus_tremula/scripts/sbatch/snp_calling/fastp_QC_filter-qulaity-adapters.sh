#!/bin/bash -l

#SBATCH -A hpc2n2025-120	
#SBATCH -n 1
#SBATCH -c 2
#SBATCH -t 1-00:00:00
#SBATCH --mail-user mimmi.eriksson@slu.se
#SBATCH --mail-type=FAIL,END
#SBATCH -J fastp
#SBATCH --output=/proj/nobackup/hpc2nstor2025-059/mimmi/aspen_snp_call/reports/sbatch/fastp/sbatch_R-%x_%j-%a.out
#SBATCH --error=/proj/nobackup/hpc2nstor2025-059/mimmi/aspen_snp_call/reports/sbatch/fastp/sbatch_R-%x_%j-%a.err
#SBATCH -a 0-39

ml GCC/12.3.0
ml GCC/13.2.0
ml fastp/0.23.4
ml multiqc/1.22.3

cd "/proj/nobackup/hpc2nstor2025-059"

data_dir="data/fastq/aspen/NordAsp/raw/reads/"
out_dir="data/fastq/aspen/NordAsp/trimmed/"
report_dir="mimmi/aspen_snp_call/reports/fastp/"

all_file_ids=($(cat mimmi/aspen_snp_call/info_files/id_fastq_files.txt))

file_id=$(echo ${all_file_ids[$SLURM_ARRAY_TASK_ID]})
sample_id=$(echo ${file_id} | cut -d"_" -f1)

fastp \
    -i ${data_dir}${file_id}_1.fq.gz \
    -I ${data_dir}${file_id}_2.fq.gz \
    -o ${out_dir}${sample_id}.R1.fq.gz \
    -O ${out_dir}${sample_id}.R2.fq.gz \
    -h ${report_dir}html/${sample_id}.html \
    -j ${report_dir}json/${sample_id}.json

