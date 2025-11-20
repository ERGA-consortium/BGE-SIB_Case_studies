#!/bin/bash

# $1 is the sample name (base name of the input file)

# Input parameters
INPUT_FILE=$1.fasta    # Input file when using genome mode
#INPUT_FILE=$1.aa       # Use this line instead for protein mode
MODE="genome"           # Analysis mode: genome, proteins, or transcriptome
LINEAGE="arthropoda_odb10"  # Lineage dataset to use
OUTPUT_NAME=$1_busco_$LINEAGE  # Output directory and file prefix
CPU=16                  # Number of CPUs to use

# Run BUSCO
busco -i "$INPUT_FILE" \
      -m "$MODE" \
      -l "$LINEAGE" \
      -o "$OUTPUT_NAME" \
      -c "$CPU" \
      --auto-lineage \
      -f

