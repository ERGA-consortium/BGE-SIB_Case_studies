#!/bin/bash

# --- Define Variables & Paths ---
ASSEMBLY="out_JBAT.FINAL.fa"       # Final assembly after manual curation
RAW_READS="all_ont_reads.fastq.gz" # Raw ONT reads for coverage
OUTPUT_PREFIX="Porphyra_nuclear"
KEEP_LIST="nuclear_keep_list.txt"
OUTPUT_FASTA="Porphyra_NUCLEAR_CLEAN.fasta"

# --- PART 1: Generating Identification Metrics ---

echo "Step 1: Running Tiara for taxonomic classification..."
tiara -i "$ASSEMBLY" -o "${OUTPUT_PREFIX}_tiara.tsv" -t 12

echo "Step 2: Calculating GC content per scaffold..."
bioawk -c fastx '{print $name, gc($seq)}' "$ASSEMBLY" > "${OUTPUT_PREFIX}_gc.tsv"

echo "Step 3: Calculating Coverage per Scaffold"
# Mapping reads to the assembly
minimap2 -ax map-ont -t 12 "$ASSEMBLY" "$RAW_READS" | samtools sort -@ 4 -o mapped.bam
samtools index mapped.bam

# General Coverage Calculation:

# # First, we calculate the Average Read Length (L) using SeqKit from the raw data statistics.
# Formula: L = Total_Bases / Total_Reads.
# Apply the general coverage formula: Coverage = (Mapped_Reads * L) / Scaffold_Length
# $3 = Mapped reads from samtools idxstats; $2 = Scaffold length.
samtools idxstats mapped.bam | \
awk -v len="$AVG_READ_LEN" 'BEGIN{OFS="\t"} {if($2>0) print $1, ($3*len)/$2}' > "${OUTPUT_PREFIX}_coverage.tsv"

# --- PART 2: Filtering and Extraction ---

# Note: You must merge the files above into a MASTER_TABLE.tsv before this step 
# or use a joined version as follows:

echo "Step 4: Filtering scaffolds based on combined metrics..."
# Logic: 
# 1. Taxonomy: Keep 'eukarya' and 'unknown' (removes bacteria/archaea)
# 2. GC Content: Keep between 0.60 (60%) and 0.85 (85%) (removes low-GC organelle scaffolds)
# 3. Coverage: Keep between 80x and 300x (removes low-cov junk and high-cov organelles/repeats)

# Assuming a MASTER_TABLE.tsv is created with columns: ID | GC | Coverage | Taxonomy
awk -F'\t' 'NR>1 && ($4 == "eukarya" || $4 == "unknown") && \
    ($2 >= 0.60 && $2 <= 0.85) && \
    ($3 >= 80 && $3 <= 300) {print $1}' MASTER_TABLE.tsv > "$KEEP_LIST"

echo "Step 5: Extracting sequences using seqtk..."
module load seqtk/1.4-GCC-13.3.0
seqtk subseq "$ASSEMBLY" "$KEEP_LIST" > "$OUTPUT_FASTA"

# Final nuclear assembly generated: $OUTPUT_FASTA
