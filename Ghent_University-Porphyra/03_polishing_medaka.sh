#!/bin/bash
# Medaka consensus polishing
# Parameters: $1=Raw Reads, $2=Racon Polished Assembly

RAW_READS=$1
INPUT_ASSEMBLY=$2
OUT_DIR="./medaka_consensus"

# Using the specified model for ONT R10.4.1 reads
medaka_consensus -i "$RAW_READS" -d "$INPUT_ASSEMBLY" -o "$OUT_DIR" -m r1041_e82_400bps_sup_v5.0.0 -t 12
