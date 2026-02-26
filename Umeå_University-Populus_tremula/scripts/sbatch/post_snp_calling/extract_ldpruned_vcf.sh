#!/bin/bash -l

#SBATCH -A hpc2n2025-120	
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 4
#SBATCH -t 1-00:00:00
#SBATCH -J extract_LD
#SBATCH --output=/home/m/mimmie/projects/aspen_hpc2nstor2025-059/mimmi/NordAsp_40/reports/sbatch/extract_LD/sbatch_R-%x_%j.out
#SBATCH --error=/home/m/mimmie/projects/aspen_hpc2nstor2025-059/mimmi/NordAsp_40/reports/sbatch/extract_LD/sbatch_R-%x_%j.err

ml GCC/13.2.0
ml BCFtools/1.19
ml HTSlib/1.19.1

cd /proj/nobackup/hpc2nstor2025-059/mimmi/NordAsp_40

input_vcf="output_data/vcf_filtered/NordAsp_biallelic_addedID.vcf.gz"
output_vcf="output_data/vcf_filtered/NordAsp_biallelic_LDpruned.vcf.gz"

# make bed of plink ld pruning output
sed 's/_/\t/g' output_data/plink/NordAsp-40_4-stairway.prune.in | \
    grep "^chr" | \
    awk '{print $1,$2-1,$2}' OFS='\t' > pruned.bed

# compress and index
bgzip pruned.bed
tabix -p bed pruned.bed.gz

# extract LD pruned snps
bcftools view \
    --threads 4 \
    --regions-file pruned.bed.gz \
    -Oz -o ${output_vcf} \
    ${input_vcf}

# index new vcf
tabix -p vcf ${output_vcf}

# clean up
rm pruned.bed.gz*