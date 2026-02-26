#!/bin/bash -l

#SBATCH -A hpc2n	
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 10
#SBATCH -t 1-00:00:00
#SBATCH -J concatVCFs
#SBATCH --output=reports/sbatch/concatVCFs/sbatch_R-%x_%j.out
#SBATCH --error=reports/sbatch/concatVCFs/sbatch_R-%x_%j.err

ml GCC/13.2.0 
ml BCFtools/1.19

bcftools concat \
    --threads 10 \
    --file-list info_files/all_vcfs_to_concat.txt \
    --write-index \
    -Oz -o data/variation/aspen/NordAsp_H204SC25041259.vcf.gz