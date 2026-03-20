#!/bin/bash -l
#SBATCH --job-name=collect_window_depths
#SBATCH --cpus-per-task=2
#SBATCH --time=04:00:00
#SBATCH --output=logs/%x-%A.slurm

ROOT="/home/people/22211215/scratch"
WD="$ROOT/tridubire"
ALIGN_DIR="$WD/02-alignment"
MOSDEPTH_DIR="$ALIGN_DIR/mosdepth"
OUT_DIR="$ALIGN_DIR/qc_summary"
mkdir -p "$OUT_DIR"

IMG_SAMTOOLS="docker://quay.io/biocontainers/samtools:1.22.1--h96c455f_0"

OUT_TSV="$OUT_DIR/window_depths.tsv"
printf "sample\tchrom\tstart\tend\tmean_depth\twindow\n" > "$OUT_TSV"

shopt -s nullglob
for f in "$MOSDEPTH_DIR"/*.win*.regions.bed.gz; do
  base=$(basename "$f")
  sample="${base%%.win*}"                              
  window=$(echo "$base" | sed -E 's/.*\.win([0-9]+).*/\1/')
  # Use samtools/htslib for zcat
  singularity exec --bind "$ROOT:$ROOT" --home "$WD" --pwd "$WD" "$IMG_SAMTOOLS" \
    bash -lc '
      zcat "'"$f"'" | awk -v s="'"$sample"'" -v w="'"$window"'" "BEGIN{OFS=\"\t\"} {print s,\$1,\$2,\$3,\$4,w}"
    ' >> "$OUT_TSV"
done

echo "Wrote: $OUT_TSV"
