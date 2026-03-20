#!/bin/bash -l

log(){ echo "[$(date +'%F %T')] $*"; }

die(){ echo "ERROR: $*" >&2; exit 1; }

run_container() {
  local image="$1"; shift
  if [[ "${CONTAINER_MODE:-enabled}" == "disabled" || -z "${image:-}" ]]; then
    "$@"
  else
    ${SINGULARITY_EXEC:-singularity exec} \
      ${SING_EXTRA_ARGS:-} \
      ${SING_BIND:+--bind "${SING_BIND}"} \
      ${SING_PWD:+--pwd "${SING_PWD}"} \
      "$image" \
      "$@"
  fi
}

require_file() {
  local f="$1"
  [[ -f "$f" ]] || die "Missing file: $f"
}

download_sra_run() {
  local run_acc="$1"
  local local_r1="${RAW_DIR}/${run_acc}_1.fastq.gz"
  local local_r2="${RAW_DIR}/${run_acc}_2.fastq.gz"

  if [[ -f "$local_r1" && -f "$local_r2" ]]; then
    log "Reads already present for $run_acc"
    return 0
  fi

  log "Downloading run $run_acc with sra-tools"

  run_container "$IMG_SRATOOLS" \
    prefetch \
      --output-directory "$SRA_CACHE_DIR" \
      "$run_acc"

  local sra_path="${SRA_CACHE_DIR}/${run_acc}/${run_acc}.sra"
  [[ -f "$sra_path" ]] || die "prefetch did not create expected file: $sra_path"

  run_container "$IMG_SRATOOLS" \
    fasterq-dump \
      --split-files \
      --threads "$THREADS" \
      --temp "$FQDUMP_TMP" \
      --outdir "$RAW_DIR" \
      "$sra_path"

  local fq1="${RAW_DIR}/${run_acc}_1.fastq"
  local fq2="${RAW_DIR}/${run_acc}_2.fastq"

  [[ -f "$fq1" ]] || die "Missing output: $fq1"
  [[ -f "$fq2" ]] || die "Missing output: $fq2"

  # Compress and remove originals
  gzip -f "$fq1"
  gzip -f "$fq2"
  rm -f "$fq1" "$fq2"

  [[ -f "$local_r1" && -f "$local_r2" ]] || \
    die "Missing gzips for $run_acc"
}

populate_samplesheet_reads_from_sra() {
  local sheet="$1"
  local tmp_sheet="${sheet}.tmp.$$"

  require_file "$sheet"

  awk -F',' '
    NR==1 {
      for (i=1; i<=NF; i++) {
        gsub(/\r/, "", $i)
        if ($i=="run_accession") run_col=i
        if ($i=="read1_filename") r1_col=i
        if ($i=="read2_filename") r2_col=i
      }
      if (!run_col || !r1_col || !r2_col) {
        print "ERROR: samplesheet must contain columns run_accession, read1_filename, read2_filename" > "/dev/stderr"
        exit 1
      }
      next
    }
    {
      gsub(/\r/, "", $run_col)
      if ($run_col != "") print $run_col
    }
  ' "$sheet" | sort -u | while read -r run_acc; do
    [[ -n "$run_acc" ]] || continue
    download_sra_run "$run_acc"
  done

  awk -F',' -v OFS=',' -v raw_dir="$RAW_DIR" '
    NR==1 {
      for (i=1; i<=NF; i++) {
        gsub(/\r/, "", $i)
        if ($i=="run_accession") run_col=i
        if ($i=="read1_filename") r1_col=i
        if ($i=="read2_filename") r2_col=i
      }
      print
      next
    }
    {
      gsub(/\r/, "", $run_col)
      $r1_col = raw_dir "/" $run_col "_1.fastq.gz"
      $r2_col = raw_dir "/" $run_col "_2.fastq.gz"
      print
    }
  ' "$sheet" > "$tmp_sheet"

  mv "$tmp_sheet" "$sheet"
  log "Update samplesheet: $sheet"
}
