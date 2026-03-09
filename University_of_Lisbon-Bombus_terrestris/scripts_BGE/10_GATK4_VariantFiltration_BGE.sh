#!/bin/bash

# Define project path
Project_path="/mnt/data3/Project_BGE_SSF_LM"

# Path to assembly (without .fna extension)
assembly_path="${Project_path}/Operational/3_reference_genomes/RepeatMasker_CHROM"
assembly="GCA_028455605.1_LMAL_v2.2.0_genomicCHROM.masked"

# Input/output directory for VCFs
vcf_dir="${Project_path}/Operational/4_data_analysis/4_snps_CHROM"
mkdir -p "$vcf_dir"

echo "🔄 Filtering SNPs for $assembly"

ref="${assembly_path}/${assembly}.fna"
raw_vcf="${vcf_dir}/cohort.${assembly}.raw.vcf.gz"
filtered_vcf="${vcf_dir}/cohort.${assembly}.flt.vcf.gz"

# Check if input VCF and reference exist
if [[ ! -f "$ref" ]]; then
  echo "❌ Reference file $ref not found. Skipping $assembly."
  exit 1
fi

if [[ ! -f "$raw_vcf" ]]; then
  echo "❌ Raw VCF $raw_vcf not found. Skipping $assembly."
  exit 1
fi

# Run VariantFiltration
gatk VariantFiltration \
  -R "$ref" \
  -V "$raw_vcf" \
  --filter-expression "QD < 2.0 || FS > 60.0 || MQ < 40.0" \
  --filter-name "basic_snp_filter" \
  -O "$filtered_vcf"

  echo "✅ Filtered VCF created: $filtered_vcf"
  echo "----------------------------------------"