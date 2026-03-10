#!/bin/bash

# Define project path
Project_path="/mnt/data3/Project_BGE_SSF_LM"

# Path to sample list
sample_list="${Project_path}/Operational/2_metadata/sample_list_Lmal.txt"

# Path to FASTQ files
fastq_path="${Project_path}/Operational/4_data_analysis/2_fastp"

# Path to assembly
assembly_path="${Project_path}/Operational/3_reference_genomes/RepeatMasker_CHROM"

# Output directories
output_base="${Project_path}/Operational/4_data_analysis/3_mapping_CHROM"
bam_dir="${output_base}/bam"
qc_dir="${output_base}/qc"

# Create output directories if they don't exist
mkdir -p "$bam_dir" "$qc_dir"

# Define single assembly (without extension)
assembly="GCA_028455605.1_LMAL_v2.2.0_genomicCHROM.masked"

# Read sample names from file
while read -r SAMPLE; do
  fastq_1="${fastq_path}/${SAMPLE}_1.fastq.gz"
  fastq_2="${fastq_path}/${SAMPLE}_2.fastq.gz"

  # Check if FASTQ files exist
  if [[ ! -f "$fastq_1" || ! -f "$fastq_2" ]]; then
    echo "❌ Missing FASTQ files for $SAMPLE. Skipping..."
    continue
  fi

  # Extract metadata from FASTQ header
  header=$(zcat "$fastq_1" | head -n 1)
  flowcell=$(echo "$header" | cut -d ':' -f3)
  lane=$(echo "$header" | cut -d ':' -f4)

  echo "🔄 Processing $SAMPLE | Flowcell: $flowcell | Lane: $lane"

  ref="${assembly_path}/${assembly}.fna"
  bam_out="${bam_dir}/${SAMPLE}.${assembly}.bam"
  qc_out="${qc_dir}/${SAMPLE}.${assembly}.flagstat.txt"

  # Check if reference file exists
  if [[ ! -f "$ref" ]]; then
    echo "❌ Reference file $ref not found. Skipping $SAMPLE."
    continue
  fi

  echo "  🧬 Aligning to $assembly..."

  # Run BWA and sort with samtools
  bwa mem -M -t 16 -R "@RG\tID:${SAMPLE}\tSM:${SAMPLE}\tPL:illumina" "$ref" "$fastq_1" "$fastq_2" | \
  samtools sort -@8 -o "$bam_out"
  samtools index "$bam_out"

  # Generate QC stats
  samtools flagstat "$bam_out" > "$qc_out"

  echo "  ✅ Done: $bam_out"
  echo "✅ Finished processing $SAMPLE"
  echo "----------------------------------------"
done < "$sample_list"
