#!/bin/bash

# Define project path
Project_path="/mnt/data3/Project_BGE_SSF_LM"

# Path to sample list
sample_list="${Project_path}/Operational/2_metadata/sample_list_Lmal.txt"

# Path to assembly (without .fna extension)
assembly_path="${Project_path}/Operational/3_reference_genomes/RepeatMasker_CHROM"
assembly="GCA_028455605.1_LMAL_v2.2.0_genomicCHROM.masked"

# Output directories
output_base="${Project_path}/Operational/4_data_analysis/3_mapping_CHROM"
bam_dir="${output_base}/bam"
qc_dir="${output_base}/qc"

# Read sample names from file
while read -r SAMPLE; do
  echo "🔄 Processing MarkDuplicates and Indexing for $SAMPLE"

  bam_in="${bam_dir}/${SAMPLE}.${assembly}.bam"
  bam_out="${bam_dir}/${SAMPLE}.${assembly}.dedup.bam"
  metrics_out="${qc_dir}/${SAMPLE}.${assembly}.dup.txt"

  # Check if input BAM exists
  if [[ ! -f "$bam_in" ]]; then
    echo "❌ BAM file $bam_in not found. Skipping $SAMPLE."
    continue
  fi

  echo "  🧬 Marking duplicates for $SAMPLE..."

  gatk MarkDuplicates \
    -I "$bam_in" \
    -O "$bam_out" \
    -M "$metrics_out"

  echo "  📌 Indexing deduplicated BAM..."

  samtools index "$bam_out"

  echo "  ✅ Done: $bam_out and index created"
  echo "----------------------------------------"
done < "$sample_list"
