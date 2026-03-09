#!/bin/bash -l

#SBATCH -A hpc2n	
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 20
#SBATCH -t 1-00:00:00
#SBATCH -J LD_stairway
#SBATCH --output=reports/sbatch/LD_stairway/sbatch_R-%x_%j.out
#SBATCH --error=reports/sbatch/LD_stairway/sbatch_R-%x_%j.err

ml GCC/10.2.0
ml PLINK/1.9b5

input_vcf="output_data/vcf_filtered/NordAsp_biallelic_addedID.vcf.gz"

# Step A – LD‑prune at r² < 0.2, 50 kb windows, sliding 10 variants
echo "LD pruning"
plink --vcf ${input_vcf} \
    --allow-extra-chr \
   	--indep-pairwise 50 10 0.2 \
  	--out output_data/plink/NordAsp-40_4-stairway