#!/bin/bash

# Define project path
Project_path="/mnt/data3/Project_BGE_SSF_LM"

# GVCF input and output directories
gvcf_dir="${Project_path}/Operational/4_data_analysis/4_snps_CHROM/gvcf"
output_dir="${Project_path}/Operational/4_data_analysis/4_snps_CHROM"
mkdir -p "$output_dir"

# Path to assembly (without .fna extension)
assembly_path="${Project_path}/Operational/3_reference_genomes/RepeatMasker_CHROM"
assembly="GCA_028455605.1_LMAL_v2.2.0_genomicCHROM.masked"

echo "🔄 Preparing list file and combining GVCFs for $assembly"

# Create list file
list_file="${gvcf_dir}/list_${assembly}.list"
> "$list_file"  # Clear existing file if any

# Find and add GVCFs to list
for gvcf in "${gvcf_dir}"/*"${assembly}.g.vcf.gz"; do
  if [[ -f "$gvcf" ]]; then
    echo "$gvcf" >> "$list_file"
  fi
done

# Check if list file has entries
if [[ ! -s "$list_file" ]]; then
  echo "❌ No GVCFs found for $assembly. Skipping CombineGVCFs."
  exit 1
fi

# Reference genome
ref="${assembly_path}/${assembly}.fna"

# Output cohort GVCF
output_gvcf="${output_dir}/cohort.${assembly}.g.vcf.gz"

# Run CombineGVCFs
gatk CombineGVCFs \
  -R "$ref" \
  -V "$list_file" \
  -O "$output_gvcf"

echo "✅ Combined GVCF created"
echo "----------------------------------------"
