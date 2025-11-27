#!/bin/bash -l

#SBATCH -A hpc2n2025-120	
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 1
#SBATCH -t 1-00:00:00
#SBATCH -J mosdepth
#SBATCH --output=/home/m/mimmie/projects/aspen_hpc2nstor2025-059/mimmi/NordAsp_40/reports/sbatch/mosdepth/sbatch_R-%x_%j-%a.out
#SBATCH --error=/home/m/mimmie/projects/aspen_hpc2nstor2025-059/mimmi/NordAsp_40/reports/sbatch/mosdepth/sbatch_R-%x_%j-%a.err
#SBATCH -a 0-39


cd /home/m/mimmie/projects/aspen_hpc2nstor2025-059/mimmi/NordAsp_40/output_data/mosdepth_coverage

eval "$(pixi shell-hook)"

bam_files=($(ls /home/m/mimmie/projects/aspen_hpc2nstor2025-059/data/populus_tremula/mapped/NordAsp/*.bam))
bam_file=$(echo ${bam_files[$SLURM_ARRAY_TASK_ID]})

sample_id=$(basename ${bam_file} _sorted_RG_dedup.bam)

mosdepth \
    --no-per-base \
    --use-median \
    --mapq 10 \
    ${sample_id} \
    ${bam_file}