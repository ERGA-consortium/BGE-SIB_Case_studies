#!/bin/bash -l
#SBATCH --job-name=bcf_sg_2n
#SBATCH --mail-user=katie.herron@ucdconnect.ie
#SBATCH --mail-type=END
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=240:00:00
#SBATCH --array=0-42%8
#SBATCH --output=logs/%x-%A_%a.slurm

# Paths
ROOT="/home/people/22211215/scratch"
WD="$ROOT/tridubire"
REF="$WD/00-data/triDubi-GCA_951804385.1.fna"

REGIONS_DIR="$WD/00-data/regions"
REG_SG1="$REGIONS_DIR/triDubi-SG1.bed"
REG_SG2="$REGIONS_DIR/triDubi-SG2.bed"

MARKDUP_DIR="$WD/02-alignment/picard_markdup"
UNIQ_DIR="$WD/02-alignment/unique_only"
OUT_DIR="$WD/03-variants/bcftools_subgen_2n"

mkdir -p "$UNIQ_DIR" "$OUT_DIR"

# Containers
IMG_SAMTOOLS="docker://quay.io/biocontainers/samtools:1.22.1--h96c455f_0"
IMG_BCFTOOLS="docker://quay.io/biocontainers/bcftools:1.22--h3a4d415_1"

# Samples list
mapfile -t SAMPLES < <(find "$MARKDUP_DIR" -maxdepth 1 -name "*.markdup.bam" \
  | sed -E 's|^.*/||; s/\.markdup\.bam$//' | sort -u)

SAMPLE="${SAMPLES[$SLURM_ARRAY_TASK_ID]}"
MARKDUP_BAM="$MARKDUP_DIR/${SAMPLE}.markdup.bam"
UNIQ_BAM="$UNIQ_DIR/${SAMPLE}.unique.bam"

OUT_SG1="$OUT_DIR/${SAMPLE}.SG1.2n.vcf.gz"
OUT_SG2="$OUT_DIR/${SAMPLE}.SG2.2n.vcf.gz"
if [[ -s "$OUT_SG1" && -s "$OUT_SG2" ]]; then
  echo "Outputs exist for $SAMPLE — skipping."
  exit 0
fi

# Ensure reference index
[[ -s "${REF}.fai" ]] || singularity exec --no-home --bind "$ROOT:$ROOT" --pwd "$WD" "$IMG_SAMTOOLS" samtools faidx "$REF"

# Step 1: create unique-only BAM if missing ---
if [[ ! -s "$UNIQ_BAM" ]]; then
  echo "Building unique-only BAM for $SAMPLE"
  singularity exec --no-home --bind "$ROOT:$ROOT" --pwd "$WD" "$IMG_SAMTOOLS" \
    bash -lc '
      # Keep: properly paired (-f 0x2); Drop: UNMAP,SECONDARY,QCFAIL,DUP,SUPPLEMENTARY (-F 0xF04); MAPQ>=40 (-q 40)
      samtools view -@ '"$SLURM_CPUS_PER_TASK"' -b -q 40 -F 0xF04 -f 0x2 "'"$MARKDUP_BAM"'" \
      | samtools sort -@ '"$SLURM_CPUS_PER_TASK"' -o "'"$UNIQ_BAM"'" -
      samtools index -@ '"$SLURM_CPUS_PER_TASK"' "'"$UNIQ_BAM"'"
    '
else
  echo "Using existing unique-only BAM: $(basename "$UNIQ_BAM")"
fi

# Diploid call by region
call_region () {
  local LABEL="$1" BED="$2"
  local OUTPFX="$OUT_DIR/${SAMPLE}.${LABEL}.2n"
  local VCF="${OUTPFX}.vcf.gz"

  [[ -s "$BED" ]] || { echo "Skip $LABEL (empty or missing $BED)"; return 0; }

  singularity exec --no-home --bind "$ROOT:$ROOT" --pwd "$WD" "$IMG_BCFTOOLS" \
    bash -lc '
      bcftools mpileup --threads '"$SLURM_CPUS_PER_TASK"' \
        -f "'"$REF"'" -R "'"$BED"'" -a AD,DP "'"$UNIQ_BAM"'" -Ou \
      | bcftools call --threads '"$SLURM_CPUS_PER_TASK"' \
        -m -Oz -o "'"$VCF"'"
      bcftools index -f "'"$VCF"'"
    '
  echo "[subgen2n] VCF: $VCF"
}

# Per-subgenome diploid calls
call_region "SG1" "$REG_SG1"
call_region "SG2" "$REG_SG2"

echo "Done for $SAMPLE"
