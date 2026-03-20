#!/bin/bash -l
#SBATCH --job-name=mosdepth_windows
#SBATCH --cpus-per-task=12
#SBATCH --time=24:00:00
#SBATCH --array=0-42%8
#SBATCH --output=logs/%x-%A_%a.slurm

# Paths
ROOT="/home/people/22211215/scratch"
WD="$ROOT/tridubire"
ALIGN_DIR="$WD/02-alignment"
MARKDUP_DIR="$ALIGN_DIR/picard_markdup"
MOSDEPTH_DIR="$ALIGN_DIR/mosdepth"
REGIONS_DIR="$WD/00-data/regions"
BED_NUC="$REGIONS_DIR/triDubi-nuclear.bed"

# Containers
IMG_MOSDEPTH="docker://quay.io/biocontainers/mosdepth:0.3.11--h0ec343a_1"
IMG_SAMTOOLS="docker://quay.io/biocontainers/samtools:1.22.1--h96c455f_0"

WINDOW="${WINDOW:-100000}"

mkdir -p "$MOSDEPTH_DIR"

# Samples
mapfile -t SAMPLES < <(find "$MARKDUP_DIR" -maxdepth 1 -name "*.markdup.bam" \
  | sed -E 's|^.*/||; s/\.markdup\.bam$//' | sort -u)

# Bounds check
[[ ${SLURM_ARRAY_TASK_ID:-0} -lt ${#SAMPLES[@]} ]] || { echo "Bad array index"; exit 1; }
SAMPLE="${SAMPLES[$SLURM_ARRAY_TASK_ID]}"

BAM="$MARKDUP_DIR/${SAMPLE}.markdup.bam"
PFX="$MOSDEPTH_DIR/${SAMPLE}.win${WINDOW}"

# Skip if already exists
if [[ -s "${PFX}.regions.bed.gz" ]]; then
  echo "Exists: $(basename "${PFX}.regions.bed.gz") — skipping."
  exit 0
fi

# A) Make windowed regions across the whole genome
singularity exec --bind "$ROOT:$ROOT" --home "$WD" --pwd "$WD" "$IMG_MOSDEPTH" \
  mosdepth --fast-mode -t "$SLURM_CPUS_PER_TASK" \
    --by "$WINDOW" \
    "$PFX" "$BAM"

# B) Filter those windows down to nuclear contigs only
NUC_LIST="${PFX}.nuc.contigs.txt"
awk '{print $1}' "$BED_NUC" | sort -u > "$NUC_LIST"

singularity exec --bind "$ROOT:$ROOT" --home "$WD" --pwd "$WD" "$IMG_SAMTOOLS" \
  bash -lc '
    zcat "'"${PFX}.regions.bed.gz"'" \
    | awk -v OFS="\t" "NR==FNR{c[\$1];next} (\$1 in c)" "'"$NUC_LIST"'" - \
    | bgzip -c > "'"${PFX}.regions.bed.gz"'.tmp" \
    && mv -f "'"${PFX}.regions.bed.gz"'.tmp" "'"${PFX}.regions.bed.gz"'" \
    && tabix -f -p bed "'"${PFX}.regions.bed.gz"'"
  '

rm -f "$NUC_LIST"

echo "Wrote: ${PFX}.regions.bed.gz   (cols: chrom  start  end  meanDepth)"
