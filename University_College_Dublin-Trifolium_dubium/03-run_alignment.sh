#!/bin/bash -l
# Align trimmed reads

ROOT="${ROOT:-/home/people/22211215/scratch}"
PROJECT_NAME="${PROJECT_NAME:-tridubire}"
WD="${WD:-${ROOT}/${PROJECT_NAME}}"

source "$WD/configs/directories.sh"
source "$WD/configs/containers.sh"
source "$WD/configs/general.sh"
source "$WD/bin/lib.sh"

FLAGSTAT_DIR="${ALIGN_DIR}/flagstat"
IDXSTAT_DIR="${ALIGN_DIR}/idxstat"

mkdir -p "$FLAGSTAT_DIR" "$IDXSTAT_DIR"

[[ -f "$SAMPLESHEET" ]] || die "Missing samplesheet: $SAMPLESHEET"
[[ -f "$REF_FASTA" ]] || die "Missing reference FASTA: $REF_FASTA"
[[ -f "${REF_FASTA}.fai" ]] || die "Missing samtools faidx index: ${REF_FASTA}.fai"
[[ -f "${REF_FASTA}.bwt.2bit.64" ]] || die "Missing bwa-mem2 index: ${REF_FASTA}.bwt.2bit.64"

col_index() {
  local needle="$1"
  awk -F',' -v needle="$needle" '
    NR==1 {
      for (i = 1; i <= NF; i++) {
        gsub(/\r/, "", $i)
        if ($i == needle) {
          print i
          exit
        }
      }
      exit 1
    }
  ' "$SAMPLESHEET"
}

RGID_COL="$(col_index "rg_id")" || die "Missing column: rg_id"
RGSM_COL="$(col_index "rg_sm")" || die "Missing column: rg_sm"
UNIT_COL="$(col_index "flowcell_lane_index")" || die "Missing column: flowcell_lane_index"
RGLB_COL="$(col_index "rg_lb")" || die "Missing column: rg_lb"
RGPL_COL="$(col_index "rg_pl")" || die "Missing column: rg_pl"
RGPU_COL="$(col_index "rg_pu")" || die "Missing column: rg_pu"

mapfile -t SAMPLES < <(
  awk -F',' -v c="$RGSM_COL" '
    NR > 1 {
      gsub(/\r/, "", $c)
      if ($c != "") print $c
    }
  ' "$SAMPLESHEET" | sort -u
)

NUM_SAMPLES=${#SAMPLES[@]}
log "Found $NUM_SAMPLES samples"

if [[ -n "${SLURM_ARRAY_TASK_ID:-}" ]]; then
  idx="${SLURM_ARRAY_TASK_ID}"
  if (( idx < 0 || idx >= NUM_SAMPLES )); then
    die "SLURM_ARRAY_TASK_ID=$idx out of range (0..$(( NUM_SAMPLES - 1 )))"
  fi
  SAMPLES=("${SAMPLES[$idx]}")
fi

for SAMPLE in "${SAMPLES[@]}"; do
  log "Running for: $SAMPLE"

  RUN_BAMS=()

  while IFS=$'\t' read -r rg_id rg_sm flowcell_lane_index rg_lb rg_pl rg_pu; do
    [[ -n "$rg_id" ]] || continue

    R1="$TRIM_DIR/${rg_sm}.${flowcell_lane_index}.R1.fq.gz"
    R2="$TRIM_DIR/${rg_sm}.${flowcell_lane_index}.R2.fq.gz"

    if [[ ! -f "$R1" || ! -f "$R2" ]]; then
      echo "WARNING: missing trimmed reads for $rg_sm / $flowcell_lane_index - skipping" >&2
      continue
    fi

    BAM_OUT="$RUN_BAM_DIR/${rg_id}.sorted.bam"
    BAI_OUT="${BAM_OUT}.bai"

    RG_STR="@RG\tID:${rg_id}\tSM:${rg_sm}\tLB:${rg_lb}\tPL:${rg_pl}\tPU:${rg_pu}"

    log "Unit: $flowcell_lane_index"
    log "  rg_id: $rg_id"
    log "  R1: $R1"
    log "  R2: $R2"

    if [[ -s "$BAM_OUT" && -s "$BAI_OUT" ]]; then
      log "Exists: $(basename "$BAM_OUT") - skipping"
    else
      log "Running BWA-MEM2 and samtools"

      run_container "$IMG_BWAMEM2" \
        bwa-mem2 mem \
          -M \
          -t "$THREADS" \
          -R "$RG_STR" \
          "$REF_FASTA" \
          "$R1" \
          "$R2" | \
      run_container "$IMG_SAMTOOLS" \
        samtools sort \
          --threads "$THREADS" \
          -o "$BAM_OUT" \
          -

      run_container "$IMG_SAMTOOLS" \
        samtools index \
          --threads "$THREADS" \
          "$BAM_OUT"
    fi

    RUN_BAMS+=("$BAM_OUT")

  done < <(
    awk -F',' \
      -v rgid_col="$RGID_COL" \
      -v rgsm_col="$RGSM_COL" \
      -v unit_col="$UNIT_COL" \
      -v rglb_col="$RGLB_COL" \
      -v rgpl_col="$RGPL_COL" \
      -v rgpu_col="$RGPU_COL" \
      -v sample="$SAMPLE" '
        NR > 1 {
          gsub(/\r/, "", $rgsm_col)
          if ($rgsm_col == sample) {
            gsub(/\r/, "", $rgid_col)
            gsub(/\r/, "", $unit_col)
            gsub(/\r/, "", $rglb_col)
            gsub(/\r/, "", $rgpl_col)
            gsub(/\r/, "", $rgpu_col)
            print $rgid_col "\t" $rgsm_col "\t" $unit_col "\t" $rglb_col "\t" $rgpl_col "\t" $rgpu_col
          }
        }
      ' "$SAMPLESHEET"
  )

  if ((${#RUN_BAMS[@]} == 0)); then
    echo "WARNING: no BAMs for $SAMPLE - skipping" >&2
    continue
  fi

  MERGED_BAM="$SAMPLE_BAM_DIR/${SAMPLE}.merged.sorted.bam"

  if [[ -s "$MERGED_BAM" && -s "${MERGED_BAM}.bai" ]]; then
    log "Merged BAM exists: $MERGED_BAM"
  else
    log "Merging ${#RUN_BAMS[@]} runs"

    if ((${#RUN_BAMS[@]} == 1)); then
      cp "${RUN_BAMS[0]}" "$MERGED_BAM"
      cp "${RUN_BAMS[0]}.bai" "${MERGED_BAM}.bai"
    else
      run_container "$IMG_SAMTOOLS" \
        samtools merge \
          --threads "$THREADS" \
          -o "$MERGED_BAM" \
          "${RUN_BAMS[@]}"

      run_container "$IMG_SAMTOOLS" \
        samtools index \
          --threads "$THREADS" \
          "$MERGED_BAM"
    fi
  fi

  log "Running samtools flagstat for $SAMPLE..."
  run_container "$IMG_SAMTOOLS" \
    samtools flagstat \
      --threads "$THREADS" \
      "$MERGED_BAM" > "$FLAGSTAT_DIR/${SAMPLE}.merged.flagstat.txt"

  log "Running samtools idxstats for $SAMPLE..."
  run_container "$IMG_SAMTOOLS" \
    samtools idxstats \
      "$MERGED_BAM" > "$IDXSTAT_DIR/${SAMPLE}.merged.idxstats.txt"

  log "Finished sample: $SAMPLE"
done

log "Alignment complete."
