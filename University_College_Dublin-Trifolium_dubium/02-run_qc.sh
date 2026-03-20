#!/bin/bash -l
# Run initial QC: FastQC, FastQ Screen, fastp trimming.

ROOT="${ROOT:-/home/people/22211215/scratch}"
PROJECT_NAME="${PROJECT_NAME:-genuspopgen}"
WD="${WD:-${ROOT}/${PROJECT_NAME}}"

source "$WD/configs/directories.sh"
source "$WD/configs/containers.sh"
source "$WD/configs/general.sh"
source "$WD/bin/lib.sh"

FASTP_HTML_DIR="${TRIM_DIR}/html"
FASTP_JSON_DIR="${TRIM_DIR}/json"

mkdir -p "$FASTP_HTML_DIR" "$FASTP_JSON_DIR"

[[ -f "$SAMPLESHEET" ]] || die "Missing samplesheet: $SAMPLESHEET"
[[ -f "$FQS_CONF" ]] || die "Missing FastQ Screen config: $FQS_CONF"

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

RGSM_COL="$(col_index "rg_sm")" || die "Missing column: rg_sm"
UNIT_COL="$(col_index "flowcell_lane_index")" || die "Missing column: flowcell_lane_index"
R1_COL="$(col_index "read1_filename")" || die "Missing column: read1_filename"
R2_COL="$(col_index "read2_filename")" || die "Missing column: read2_filename"

mapfile -t SAMPLES < <(
  awk -F',' -v c="$RGSM_COL" '
    NR > 1 {
      gsub(/\r/, "", $c)
      if ($c != "") print $c
    }
  ' "$SAMPLESHEET" | sort -u
)

if [[ -n "${SLURM_ARRAY_TASK_ID:-}" ]]; then
  idx="${SLURM_ARRAY_TASK_ID}"
  if (( idx < 0 || idx >= ${#SAMPLES[@]} )); then
    die "SLURM_ARRAY_TASK_ID=$idx out of range (0..$(( ${#SAMPLES[@]} - 1 )))"
  fi
  SAMPLES=("${SAMPLES[$idx]}")
fi

log "Initial QC for ${#SAMPLES[@]} sample(s)."
log "Using samplesheet: $SAMPLESHEET"
log "Using FastQ Screen config: $FQS_CONF"
log "FastQC output dir: $RAW_FASTQC_DIR"
log "FastQ Screen output dir: $FQS_OUT_DIR"
log "Trimmed reads dir: $TRIM_DIR"

for SAMPLE in "${SAMPLES[@]}"; do
  log "Sample: $SAMPLE"

  awk -F',' \
    -v rgsm_col="$RGSM_COL" \
    -v unit_col="$UNIT_COL" \
    -v r1_col="$R1_COL" \
    -v r2_col="$R2_COL" \
    -v sample="$SAMPLE" '
      NR > 1 {
        gsub(/\r/, "", $rgsm_col)
        if ($rgsm_col == sample) {
          gsub(/\r/, "", $unit_col)
          gsub(/\r/, "", $r1_col)
          gsub(/\r/, "", $r2_col)
          print $unit_col "\t" $r1_col "\t" $r2_col
        }
      }
    ' "$SAMPLESHEET" |
  while IFS=$'\t' read -r unit r1 r2; do
    [[ -n "$r1" ]] || die "Blank read1_filename for sample $SAMPLE unit $unit"
    [[ -n "$r2" ]] || die "Blank read2_filename for sample $SAMPLE unit $unit"

    [[ -f "$r1" ]] || die "Missing read1 file: $r1"
    [[ -f "$r2" ]] || die "Missing read2 file: $r2"

    trim_r1="${TRIM_DIR}/${SAMPLE}.${unit}.R1.fq.gz"
    trim_r2="${TRIM_DIR}/${SAMPLE}.${unit}.R2.fq.gz"
    fastp_html="${FASTP_HTML_DIR}/${SAMPLE}.${unit}.html"
    fastp_json="${FASTP_JSON_DIR}/${SAMPLE}.${unit}.json"

    log "Unit: $unit"
    log "  R1: $r1"
    log "  R2: $r2"

    # FastQC on raw reads
    run_container "$IMG_FASTQC" fastqc \
      -t "$THREADS" \
      -o "$RAW_FASTQC_DIR" \
      "$r1" "$r2"

    # FastQ Screen on raw reads
    run_container "$IMG_FASTQ_SCREEN" fastq_screen \
      --aligner bwa \
      --threads "$THREADS" \
      --subset 750000 \
      --conf "$FQS_CONF" \
      --outdir "$FQS_OUT_DIR" \
      "$r1" "$r2"

    # Trim reads per unit
    if [[ -f "$trim_r1" && -f "$trim_r2" ]]; then
      log "Skipping fastp - trimmed reads already exist."
    else
      log "Running fastp..."
      run_container "$IMG_FASTP" fastp \
        -w "$THREADS" \
        -i "$r1" \
        -I "$r2" \
        -o "$trim_r1" \
        -O "$trim_r2" \
        -h "$fastp_html" \
        -j "$fastp_json"
    fi
  done
done

log "Initial QC + trimming complete."
