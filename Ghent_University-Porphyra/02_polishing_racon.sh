#!/bin/bash
# 3-round Racon polishing loop
# Parameters: $1=Raw Reads, $2=Initial Assembly

RAW_READS=$1
CURRENT_ASSEMBLY=$2
OUT_DIR="./racon_polished"
mkdir -p $OUT_DIR

for i in {1..3}
do
    echo "Starting Racon Round $i"
    # Mapping
    minimap2 -x map-ont -t 16 $CURRENT_ASSEMBLY $RAW_READS > overlaps.paf
    # Polishing
    racon -t 16 $RAW_READS overlaps.paf $CURRENT_ASSEMBLY > "${OUT_DIR}/polished_v${i}.fasta"
    # Update for next round
    CURRENT_ASSEMBLY="${OUT_DIR}/polished_v${i}.fasta"
    rm overlaps.paf
done
