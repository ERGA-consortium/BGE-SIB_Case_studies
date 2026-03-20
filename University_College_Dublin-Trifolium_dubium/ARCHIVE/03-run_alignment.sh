#!/bin/bash -l
#SBATCH --job-name=alignment
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
TRIM_DIR="$WD/01-qc/fastp_trimmed"
ALIGN_DIR="$WD/02-alignment"
RUN_BAM_DIR="$ALIGN_DIR/bam-per_run"
SAMPLE_BAM_DIR="$ALIGN_DIR/bam-merged"
FLAGSTAT_DIR="$ALIGN_DIR/flagstat"
IDXSTAT_DIR="$ALIGN_DIR/idxstat"
REF_GENOME="$ROOT/tridubire/00-data/triDubi-GCA_951804385.1.fna"

mkdir -p "$RUN_BAM_DIR" "$SAMPLE_BAM_DIR" "$FLAGSTAT_DIR" "$IDXSTAT_DIR"

# Singularity images
IMG_BWAMEM2="docker://quay.io/biocontainers/bwa-mem2:2.3--he70b90d_0"
IMG_SAMTOOLS="docker://quay.io/biocontainers/samtools:1.22.1--h96c455f_0"

# Get sample
mapfile -t SAMPLES < <(
  find "$TRIM_DIR" -maxdepth 1 -name "*_R1.trimmed.fq.gz" \
  | sed -E 's|^.*/||' \
  | awk -F'_' '{print $1}' \
  | sort -u
)

NUM_SAMPLES=${#SAMPLES[@]}
echo "Found $NUM_SAMPLES unique samples."

SAMPLE="${SAMPLES[$SLURM_ARRAY_TASK_ID]}"
echo "Running alignment for sample: $SAMPLE"

mapfile -t RUNS < <(
  find "$TRIM_DIR" -maxdepth 1 -name "${SAMPLE}_*_R1.trimmed.fq.gz" \
  | sed -E 's|^.*/||' \
  | sed -E 's/_R1\.trimmed\.fq\.gz$//' \
  | awk -F'_' '{print $2"."$4"."$3}' \
  | sort -u
)

NUM_RUNS=${#RUNS[@]}
echo "Found $NUM_RUNS runs for $SAMPLE."

RUN_BAMS=()

for RUN in "${RUNS[@]}"; do
  IFS='.' read -r FLOWCELL LANE INDEX <<< "$RUN"
  RUN_TAG="${FLOWCELL}.${LANE}.${INDEX}"
  echo "Run: $RUN_TAG"

  R1="$TRIM_DIR/${SAMPLE}_${FLOWCELL}_${INDEX}_${LANE}_R1.trimmed.fq.gz"
  R2="$TRIM_DIR/${SAMPLE}_${FLOWCELL}_${INDEX}_${LANE}_R2.trimmed.fq.gz"

  RG_STR="@RG\tID:${SAMPLE}.${FLOWCELL}.${LANE}.${INDEX}\tSM:${SAMPLE}\tLB:${SAMPLE}_lib1\tPL:ILLUMINA\tPU:${FLOWCELL}.${LANE}.${INDEX}"

  SAM_OUT="$RUN_BAM_DIR/${SAMPLE}.${FLOWCELL}.${LANE}.${INDEX}.sam"
  BAM_OUT="$RUN_BAM_DIR/${SAMPLE}.${FLOWCELL}.${LANE}.${INDEX}.sorted.bam"

  if [[ -s "$BAM_OUT" ]]; then
    echo "Exists: $(basename "$BAM_OUT") - skipping."
  else
    echo "Running BWA-MEM2 for $SAMPLE"
    singularity exec --no-home --bind "$ROOT:$ROOT" --pwd "$WD" "$IMG_BWAMEM2" \
      bwa-mem2 mem -M -t "$SLURM_CPUS_PER_TASK" -R "$RG_STR" "$REF_GENOME" "$R1" "$R2" \
      > "$SAM_OUT"

    echo "Running samtools sort for $SAMPLE"
    singularity exec --no-home --bind "$ROOT:$ROOT" --pwd "$WD" "$IMG_SAMTOOLS" \
      samtools sort --threads "$SLURM_CPUS_PER_TASK" -o "$BAM_OUT" "$SAM_OUT"
  fi
  RUN_BAMS+=("$BAM_OUT")
done

MERGED_BAM="$SAMPLE_BAM_DIR/${SAMPLE}.merged.sorted.bam"

if [[ -s "$MERGED_BAM" ]]; then
  echo "Merged BAM exists: $MERGED_BAM"
else
  echo "Merging ${#RUN_BAMS[@]} runs"
  singularity exec --no-home --bind "$ROOT:$ROOT" --pwd "$WD" "$IMG_SAMTOOLS" \
    samtools merge --threads "$SLURM_CPUS_PER_TASK" -o "$MERGED_BAM" "${RUN_BAMS[@]}"
  singularity exec --no-home --bind "$ROOT:$ROOT" --pwd "$WD" "$IMG_SAMTOOLS" \
    samtools index "$MERGED_BAM"
fi

echo "Running samtools flagstat for $SAMPLE"
singularity exec --no-home --bind "$ROOT:$ROOT" --pwd "$WD" "$IMG_SAMTOOLS" \
  samtools flagstat --threads "$SLURM_CPUS_PER_TASK" "$MERGED_BAM" > "$FLAGSTAT_DIR/${SAMPLE}.merged.flagstat.txt"

echo "Running samtools idxstat for $SAMPLE"
singularity exec --no-home --bind "$ROOT:$ROOT" --pwd "$WD" "$IMG_SAMTOOLS" \
  samtools idxstats "$MERGED_BAM" > "$IDXSTAT_DIR/${SAMPLE}.merged.idxstats.txt"

echo "Finished sample: $SAMPLE"
