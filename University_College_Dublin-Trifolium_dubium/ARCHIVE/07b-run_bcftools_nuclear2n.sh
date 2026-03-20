#!/bin/bash -l
#SBATCH --job-name=bcf_nuclear_2n
#SBATCH --cpus-per-task=16
#SBATCH --time=48:00:00
#SBATCH --array=0-42%8
#SBATCH --output=logs/%x-%A_%a.slurm

# Paths
ROOT="/home/people/22211215/scratch"
WD="$ROOT/tridubire"
REF="$WD/00-data/triDubi-GCA_951804385.1.fna"

REGIONS_DIR="$WD/00-data/regions"
REG_NUC="$REGIONS_DIR/triDubi-nuclear.bed"

MARKDUP_DIR="$WD/02-alignment/picard_markdup"
UNIQ_DIR="$WD/02-alignment/unique_only"
OUT_DIR="$WD/03-variants/bcftools_nuclear_2n"

mkdir -p "$UNIQ_DIR" "$OUT_DIR"

# Containers
IMG_SAMTOOLS="docker://quay.io/biocontainers/samtools:1.22.1--h96c455f_0"
IMG_BCFTOOLS="docker://quay.io/biocontainers/bcftools:1.22--h3a4d415_1"

# Samples
mapfile -t SAMPLES < <(find "$MARKDUP_DIR" -maxdepth 1 -name "*.markdup.bam" \
  | sed -E 's|^.*/||; s/\.markdup\.bam$//' | sort -u)

SAMPLE="${SAMPLES[$SLURM_ARRAY_TASK_ID]}"

MARKDUP_BAM="$MARKDUP_DIR/${SAMPLE}.markdup.bam"
UNIQ_BAM="$UNIQ_DIR/${SAMPLE}.unique.bam"
VCF="$OUT_DIR/${SAMPLE}.nuclear.2n.vcf.gz"

OUT_NUC="$OUT_DIR/${SAMPLE}.nuclear.2n.vcf.gz"
if [[ -s "$OUT_NUC" ]]; then
  echo "Output exists for $SAMPLE — skipping."
  exit 0
fi

[[ -s "${REF}.fai" ]] || singularity exec --no-home --bind "$ROOT:$ROOT" --pwd "$WD" "$IMG_SAMTOOLS" samtools faidx "$REF"

if [[ ! -s "$UNIQ_BAM" ]]; then
  echo "Building unique-only BAM for $SAMPLE"
  singularity exec --no-home --bind "$ROOT:$ROOT" --pwd "$WD" "$IMG_SAMTOOLS" \
    bash -lc '
      samtools view -@ '"$SLURM_CPUS_PER_TASK"' -b -q 40 -F 0xF04 -f 0x2 "'"$MARKDUP_BAM"'" \
      | samtools sort -@ '"$SLURM_CPUS_PER_TASK"' -o "'"$UNIQ_BAM"'" -
      samtools index -@ '"$SLURM_CPUS_PER_TASK"' "'"$UNIQ_BAM"'"
    '
else
  echo "Using existing unique-only BAM: $(basename "$UNIQ_BAM")"
fi

singularity exec --no-home --bind "$ROOT:$ROOT" --pwd "$WD" "$IMG_BCFTOOLS" \
  bash -lc '
    bcftools mpileup --threads '"$SLURM_CPUS_PER_TASK"' \
      -f "'"$REF"'" -R "'"$REG_NUC"'" -a AD,DP "'"$UNIQ_BAM"'" -Ou \
    | bcftools call --threads '"$SLURM_CPUS_PER_TASK"' \
      -m -Oz -o "'"$VCF"'"
    bcftools index -f "'"$VCF"'"
  '

echo "[nuclear2n] VCF: $VCF"
