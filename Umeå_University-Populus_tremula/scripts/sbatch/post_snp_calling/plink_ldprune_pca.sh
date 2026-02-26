#!/bin/bash -l

#SBATCH -A hpc2n2025-120	
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 1
#SBATCH -t 1-00:00:00
#SBATCH -J ldprune_pca
#SBATCH --output=/home/m/mimmie/projects/aspen_hpc2nstor2025-059/mimmi/NordAsp_40/reports/sbatch/ldprune_pca/sbatch_R-%x_%j.out
#SBATCH --error=/home/m/mimmie/projects/aspen_hpc2nstor2025-059/mimmi/NordAsp_40/reports/sbatch/ldprune_pca/sbatch_R-%x_%j.err

ml GCC/10.2.0
ml PLINK/1.9b5

cd /proj/nobackup/hpc2nstor2025-059/mimmi/NordAsp_40

input_vcf="output_data/vcf_filtered/NordAsp_H204SC25041259_biallelic-missing0.8-maf0.05_recoded.vcf.gz"

# Step A – LD‑prune at r² < 0.2, 50 kb windows, sliding 10 variants
echo "LD pruning"
plink --vcf ${input_vcf} \
    --allow-extra-chr \
   	--indep-pairwise 50 10 0.2 \
  	--out output_data/plink/NordAsp_40

# Step B make freq file
echo "create a .freq file"
plink --vcf ${input_vcf} \
    --allow-extra-chr \
    --extract output_data/plink/NordAsp_40.prune.in \
    --freq \
  	--out output_data/plink/NordAsp_40

# Principal component analysis 
# Step C - PCA on the pruned SNP subset
echo "run plink pca"
plink --vcf ${input_vcf} \
    --allow-extra-chr \
    --extract output_data/plink/NordAsp_40.prune.in \
    --read-freq output_data/plink/NordAsp_40.frq \
   	--pca \
   	--out output_data/plink/NordAsp_40_pca
