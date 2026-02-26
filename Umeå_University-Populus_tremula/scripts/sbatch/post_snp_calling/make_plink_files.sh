#!/bin/bash -l

#SBATCH -A hpc2n	
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 1
#SBATCH -t 1-00:00:00
#SBATCH -J make_plink
#SBATCH --output=reports/sbatch/make_plink/sbatch_R-%x_%j.out
#SBATCH --error=reports/sbatch/make_plink/sbatch_R-%x_%j.err

ml GCC/10.2.0
ml PLINK/1.9b5

input_vcf="output_data/vcf_filtered/NordAsp_H204SC25041259_biallelic-missing0.8-maf0.05_recoded.vcf.gz"

# make plink files for filtered vcf
plink \
    --make-bed \
    --allow-extra-chr \
    --vcf ${input_vcf} \
    --out output_data/plink/NordAsp_40

# make plink files with the LD pruned set
plink \
    --make-bed \
    --allow-extra-chr \
    --vcf ${input_vcf} \
    --extract output_data/plink/NordAsp_40.prune.in \
    --out output_data/plink/NordAsp_40_LDpruned