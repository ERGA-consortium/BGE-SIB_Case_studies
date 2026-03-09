#!/bin/bash -l

#SBATCH -A hpc2n	
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 1
#SBATCH -t 1-00:00:00
#SBATCH -J nt_div
#SBATCH --output=reports/sbatch/nt_div/sbatch_R-%x_%j.out
#SBATCH --error=reports/sbatch/nt_div/sbatch_R-%x_%j.err

ml GCC/13.2.0
ml VCFtools/0.1.16

input_vcf="output_data/vcf_filtered/NordAsp_H204SC25041259_biallelic-missing0.8-maf0.05_recoded.vcf.gz"

vcftools --gzvcf ${input_vcf} \
     	--site-pi \
     	--out output_data/vcftools/NordAsp40_pi