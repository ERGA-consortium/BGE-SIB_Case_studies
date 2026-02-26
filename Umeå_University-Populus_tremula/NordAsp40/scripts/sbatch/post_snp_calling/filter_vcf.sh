#!/bin/bash -l

#SBATCH -A hpc2n2025-120	
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 1
#SBATCH -t 1-00:00:00
#SBATCH -J filter_vcf
#SBATCH --output=/home/m/mimmie/projects/aspen_hpc2nstor2025-059/mimmi/NordAsp_40/reports/sbatch/filter_vcf/sbatch_R-%x_%j.out
#SBATCH --error=/home/m/mimmie/projects/aspen_hpc2nstor2025-059/mimmi/NordAsp_40/reports/sbatch/filter_vcf/sbatch_R-%x_%j.err

ml GCC/13.2.0
ml BCFtools/1.19
ml VCFtools/0.1.16
ml HTSlib/1.19.1

cd /proj/nobackup/hpc2nstor2025-059/mimmi/NordAsp_40

input_vcf="/proj/nobackup/hpc2nstor2025-059/data/populus_tremula/variation/NordAsp/NordAsp_H204SC25041259.vcf.gz"
# output_vcf="output_data/vcf_filtered/NordAsp_biallelic-missing0.8-maf0.05_recoded.vcf.gz"
output_vcf="output_data/vcf_filtered/NordAsp_biallelic_addedID.vcf.gz"

# filter vcf file
bcftools view \
    -m2 \
    -M2 \
    -v snps \
    ${input_vcf} | \
vcftools --vcf - \
    --max-missing 0.8 \
    --maf 0.05 \
    --min-alleles 2 \
    --max-alleles 2 \
    --recode \
    --recode-INFO-all \
    --stdout | \
bcftools view -Oz \
    -o output_data/vcf_filtered/temp.vcf.gz

# index vcf
tabix -p vcf output_data/vcf_filtered/temp.vcf.gz

# modify ID column to be "chr1_pos_ref_alt"
# produce a bed file that contains the new ID format
bcftools view --no-version -H output_data/vcf_filtered/temp.vcf.gz | \
awk '{print $1,$2-1,$2,$1"_"$2"_"$4"_"$5}' OFS='\t' > annot_ID.bed

# compress it and index it
bgzip annot_ID.bed
tabix -p bed annot_ID.bed.gz

# annotate ID column with the new format
bcftools annotate \
    -c CHROM,FROM,TO,ID \
    -a annot_ID.bed.gz \
    -Oz -o ${output_vcf} \
    output_data/vcf_filtered/temp.vcf.gz

tabix -p vcf ${output_vcf}

# clean up
rm output_data/vcf_filtered/temp.vcf.gz*
rm annot_ID.bed.gz*