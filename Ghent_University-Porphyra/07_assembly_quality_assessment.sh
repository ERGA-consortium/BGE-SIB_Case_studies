#!/bin/bash

# --- Define Variables ---
ASSEMBLY="Porphyra_NUCLEAR_CLEAN.fasta"
THREADS=12
OUT_DIR_QUAST="./quast_results"
OUT_DIR_BUSCO="./busco_results"

# 1. QUAST Analysis (Continuity Metrics)
# QUAST calculates N50, L50, total length, and GC content.
echo "Running QUAST for assembly continuity assessment..."
quast.py "$ASSEMBLY" -o "$OUT_DIR_QUAST" -t "$THREADS"

# 2. BUSCO Analysis (Biological Completeness)
# We run BUSCO against two databases to ensure thorough validation:
# - Eukaryota: Broad evolutionarily conserved genes.
# - Rhodophyta: Specific genes conserved across Red Algae.

# Run BUSCO with Eukaryota database (eukaryota_odb10)
echo "Running BUSCO with Eukaryota database..."
busco -i "$ASSEMBLY" \
      -l eukaryota_odb10 \
      -o busco_euk \
      -m genome --cpu "$THREADS" \
      --out_path "$OUT_DIR_BUSCO"

# Run BUSCO with Rhodophyta database (rhodophyta_odb12)
echo "Running BUSCO with Rhodophyta database..."
busco -i "$ASSEMBLY" \
      -l rhodophyta_odb12 \
      -o busco_rhodo \
      -m genome --cpu "$THREADS" \
      --out_path "$OUT_DIR_BUSCO"

echo "Quality Assessment complete. Results are in $OUT_DIR_QUAST and $OUT_DIR_BUSCO"
