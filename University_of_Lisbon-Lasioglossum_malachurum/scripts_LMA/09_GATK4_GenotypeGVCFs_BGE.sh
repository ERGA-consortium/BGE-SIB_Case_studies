#!/bin/bash

# Define project path
Project_path="/mnt/data3/Project_BGE_SSF_LM"

# Path to assembly (without .fna extension)
assembly_path="${Project_path}/Operational/3_reference_genomes/RepeatMasker_CHROM"
assembly="GCA_028455605.1_LMAL_v2.2.0_genomicCHROM.masked"

# Input/output directory for cohort GVCFs and VCFs
cohort_dir="${Project_path}/Operational/4_data_analysis/4_snps_CHROM"
mkdir -p "$cohort_dir"

echo "🔄 Genotyping cohort for $assembly"

ref="${assembly_path}/${assembly}.fna"
cohort_gvcf="${cohort_dir}/cohort.${assembly}.g.vcf.gz"
output_vcf="${cohort_dir}/cohort.${assembly}.raw.vcf.gz"

# Check if input GVCF and reference exist
if [[ ! -f "$ref" ]]; then
  echo "❌ Reference file $ref not found. Skipping $assembly."
  exit 1
fi

if [[ ! -f "$cohort_gvcf" ]]; then
  echo "❌ Cohort GVCF $cohort_gvcf not found. Skipping $assembly."
  exit 1
fi

# Run GenotypeGVCFs
gatk GenotypeGVCFs \
  -R "$ref" \
  -V "$cohort_gvcf" \
  -O "$output_vcf"

echo "✅ Raw VCF created: $output_vcf"
echo "----------------------------------------"