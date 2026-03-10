#!/bin/bash

# Define project path
Project_path="/mnt/data3/Project_BGE_SSF_LM"

# Path to sample list
sample_list="${Project_path}/Operational/2_metadata/sample_list_Lmal.txt"

# Path to assembly (without .fna extension)
assembly_path="${Project_path}/Operational/3_reference_genomes/RepeatMasker_CHROM"
assembly="GCA_028455605.1_LMAL_v2.2.0_genomicCHROM.masked"

# Path to deduplicated BAM files
bam_dir="${Project_path}/Operational/4_data_analysis/3_mapping_CHROM/bam"

# Output directory for GVCFs
gvcf_dir="${Project_path}/Operational/4_data_analysis/4_snps_CHROM/gvcf"
mkdir -p "$gvcf_dir"

# Read sample names from file
while read -r SAMPLE; do
  echo "🔄 Running HaplotypeCaller for $SAMPLE"

  ref="${assembly_path}/${assembly}.fna"
  bam_in="${bam_dir}/${SAMPLE}.${assembly}.dedup.bam"
  gvcf_out="${gvcf_dir}/${SAMPLE}.${assembly}.g.vcf.gz"

  # Check if input BAM and reference exist
  if [[ ! -f "$bam_in" ]]; then
    echo "❌ Deduplicated BAM file $bam_in not found. Skipping $SAMPLE."
    continue
  fi

  if [[ ! -f "$ref" ]]; then
    echo "❌ Reference file $ref not found. Skipping $SAMPLE."
    continue
  fi

  echo "  🧬 Calling variants for $SAMPLE..."

  gatk HaplotypeCaller \
    -R "$ref" \
    -I "$bam_in" \
    -ERC GVCF \
    -O "$gvcf_out"

  echo "  ✅ GVCF created: $gvcf_out"
  echo "----------------------------------------"
done < "$sample_list"
