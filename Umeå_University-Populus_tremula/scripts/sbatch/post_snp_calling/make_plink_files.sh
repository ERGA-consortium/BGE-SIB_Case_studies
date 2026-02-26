#!/bin/bash -l

#SBATCH -A hpc2n2025-120	
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 1
#SBATCH -t 1-00:00:00
#SBATCH -J make_plink
#SBATCH --output=/home/m/mimmie/projects/aspen_hpc2nstor2025-059/mimmi/NordAsp_40/reports/sbatch/make_plink/sbatch_R-%x_%j.out
#SBATCH --error=/home/m/mimmie/projects/aspen_hpc2nstor2025-059/mimmi/NordAsp_40/reports/sbatch/make_plink/sbatch_R-%x_%j.err

ml GCC/10.2.0
ml PLINK/1.9b5

cd /proj/nobackup/hpc2nstor2025-059/mimmi/NordAsp_40

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