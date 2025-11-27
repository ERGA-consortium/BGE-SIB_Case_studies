#!/bin/bash -l

#SBATCH -A hpc2n2025-120	
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 20
#SBATCH -t 1-00:00:00
#SBATCH -J LD_stairway
#SBATCH --output=/home/m/mimmie/projects/aspen_hpc2nstor2025-059/mimmi/NordAsp_40/reports/sbatch/LD_stairway/sbatch_R-%x_%j.out
#SBATCH --error=/home/m/mimmie/projects/aspen_hpc2nstor2025-059/mimmi/NordAsp_40/reports/sbatch/LD_stairway/sbatch_R-%x_%j.err

ml GCC/10.2.0
ml PLINK/1.9b5

cd /proj/nobackup/hpc2nstor2025-059/mimmi/NordAsp_40

input_vcf="output_data/vcf_filtered/NordAsp_biallelic_addedID.vcf.gz"

# Step A – LD‑prune at r² < 0.2, 50 kb windows, sliding 10 variants
echo "LD pruning"
plink --vcf ${input_vcf} \
    --allow-extra-chr \
   	--indep-pairwise 50 10 0.2 \
  	--out output_data/plink/NordAsp-40_4-stairway