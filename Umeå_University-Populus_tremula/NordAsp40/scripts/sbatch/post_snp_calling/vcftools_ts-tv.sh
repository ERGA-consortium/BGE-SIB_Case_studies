#!/bin/bash -l

#SBATCH -A hpc2n2025-120	
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 1
#SBATCH -t 1-00:00:00
#SBATCH -J vcftools_ts-tv
#SBATCH --output=/home/m/mimmie/projects/aspen_hpc2nstor2025-059/mimmi/NordAsp_40/reports/sbatch/vcftools_ts-tv/sbatch_R-%x_%j.out
#SBATCH --error=/home/m/mimmie/projects/aspen_hpc2nstor2025-059/mimmi/NordAsp_40/reports/sbatch/vcftools_ts-tv/sbatch_R-%x_%j.err

ml GCC/13.2.0
ml VCFtools/0.1.16

cd /proj/nobackup/hpc2nstor2025-059/mimmi/NordAsp_40

# Summarise transition transverion ratio
vcftools --gzvcf output_data/vcf_filtered/NordAsp_biallelic_addedID.vcf.gz \
    --TsTv-summary \
    --out output_data/vcftools/NordAsp40
