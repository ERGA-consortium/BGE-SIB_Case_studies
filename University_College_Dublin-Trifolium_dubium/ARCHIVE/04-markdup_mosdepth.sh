#!/bin/bash -l
#SBATCH --job-name=markdup_mosdepth
#SBATCH --mail-user=katie.herron@ucdconnect.ie
#SBATCH --mail-type=END
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=240:00:00
#SBATCH --array=0-42%8
#SBATCH --output=logs/%x-%A_%a.slurm

# Directories
ROOT="/home/people/22211215/scratch/"
WD="$ROOT/tridubire"
ALIGN_DIR="$ROOT/tridubire/02-alignment"
SAMPLE_BAM_DIR="$ALIGN_DIR/bam-merged"
MARKDUP_DIR="$ALIGN_DIR/picard_markdup"
FLAGSTAT_DIR="$ALIGN_DIR/flagstat"
IDXSTAT_DIR="$ALIGN_DIR/idxstat"
MOSDEPTH_DIR="$ALIGN_DIR/mosdepth"

# Singularity images
IMG_SAMTOOLS="docker://quay.io/biocontainers/samtools:1.22.1--h96c455f_0"
IMG_PICARD="docker://quay.io/biocontainers/picard:3.4.0--hdfd78af_0"
IMG_MOSDEPTH="docker://quay.io/biocontainers/mosdepth:0.3.11--h0ec343a_1"

mkdir -p "$ALIGN_DIR" "$MARKDUP_DIR" "$FLAGSTAT_DIR" "$IDXSTAT_DIR" "$MOSDEPTH_DIR"

# Get sample
mapfile -t SAMPLES < <(
  find "$SAMPLE_BAM_DIR" -maxdepth 1 -name "*.merged.sorted.bam" \
  | sed -E 's|^.*/||' \
  | sed -E 's/\.merged\.sorted\.bam$//' \
  | sort -u
)

NUM_SAMPLES=${#SAMPLES[@]}
echo "Found $NUM_SAMPLES unique samples."

SAMPLE="${SAMPLES[$SLURM_ARRAY_TASK_ID]}"
echo "Marking duplicates for sample: $SAMPLE"

# 1. Mark duplicates
singularity exec --no-home --bind "$ROOT:$ROOT" --pwd "$WD" "$IMG_PICARD" \
  picard MarkDuplicates \
  I="$SAMPLE_BAM_DIR/${SAMPLE}.merged.sorted.bam" \
  O="$MARKDUP_DIR/${SAMPLE}.markdup.bam" \
  M="$MARKDUP_DIR/${SAMPLE}.dup_metrics.txt" \
  VALIDATION_STRINGENCY=LENIENT CREATE_INDEX=true

# 2. Flagstat post-duplicate marking
singularity exec --no-home --bind "$ROOT:$ROOT" --pwd "$WD" "$IMG_SAMTOOLS" \
  samtools flagstat --threads "$SLURM_CPUS_PER_TASK" "$MARKDUP_DIR/${SAMPLE}.markdup.bam" \
  > "$FLAGSTAT_DIR/${SAMPLE}.markdup.flagstat.txt"

# 3. Idxstats post-duplicate marking
singularity exec --no-home --bind "$ROOT:$ROOT" --pwd "$WD" "$IMG_SAMTOOLS" \
  samtools idxstats "$MARKDUP_DIR/${SAMPLE}.markdup.bam" \
  > "$IDXSTAT_DIR/${SAMPLE}.markdup.idxstats.txt"

# 4. Per-sample depth with mosdepth
singularity exec --no-home --bind "$ROOT:$ROOT" --pwd "$WD" "$IMG_MOSDEPTH" \
  mosdepth --fast-mode --threads "$SLURM_CPUS_PER_TASK" \
    -n "$MOSDEPTH_DIR/${SAMPLE}" "$MARKDUP_DIR/${SAMPLE}.markdup.bam"

echo "Finished sample: $SAMPLE"
