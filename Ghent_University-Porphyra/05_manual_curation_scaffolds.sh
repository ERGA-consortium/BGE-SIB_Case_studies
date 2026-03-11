#!/bin/bash

# 1. Explanation of Juicebox Curation
# After automated scaffolding with HapHiC, the .hic and .assembly files 
# were loaded into Juicebox (v1.8.8) to:
# - Correct misjoins and translocations.
# - Reorient scaffolds based on Hi-C contact signal.
# - Assign final chromosome boundaries.

# --- Define Variables ---
JUICER_POST="/home/jmorcill/HapHiC/utils/juicer"
ASSEMBLY_FA="consensus.fasta" # polished input assembly
REVIEW_FILE="out_JBAT.review.assembly" # Generated after Juicebox curation
LIFTOVER_AGP="out_JBAT.liftover.agp"

# 2. Generate the final FASTA file for the scaffolds
# This command applies the manual changes to the actual sequence.
echo "Generating final curated FASTA file..."
$JUICER_POST post \
    -o final_curated_assembly \
    "$REVIEW_FILE" \
    "$LIFTOVER_AGP" \
    "$ASSEMBLY_FA"

# Final chromosome-level assembly: final_curated_assembly.fasta"
