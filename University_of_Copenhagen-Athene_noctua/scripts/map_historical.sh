#!/usr/bin/env bash
#SBATCH --job-name=historical_bwa
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.err
#SBATCH --cpus-per-task=4
#SBATCH --mem=20G
#SBATCH --time=00:20:00

set -euo pipefail
# Report both sample & lane on any error
trap 'echo "[$(date "+%F %T")] ❌ Error during processing ${SAMPLE}_${LANE}" >&2' ERR

timestamp(){ date '+%F %T'; }

# Load modules
module load bwa samtools mapdamage2/2.2.2 gcc R/3.6.1 \
            fastp/0.24.0 perl openjdk fastqc/0.12.1

# SeqPrep2 path
SEQPREP=/projects/erode/people/nrv690/nrv690.hernan/software/SeqPrep2/SeqPrep2

if [[ $# -ne 6 ]]; then
  echo "Usage: $0 <sample_id> <manifest.tsv> <reference.fasta> <ref_name> <threads> <outdir>"
  exit 1
fi

SAMPLE=$1
MANIFEST=$2
REFERENCE=$3
REF_NAME=$4
THREADS=$5
OUTDIR=$6

MAX_MISINCORP_FREQUENCY=0.4
READ_PLOT_LENGTH=50
MAPDOWNSAMPLE=10000

mkdir -p "$OUTDIR" tmp trimmed out_fastqc out_fastqc_trimmed out_mapdamage
exec > >(tee -i "${OUTDIR}/pipeline_${SAMPLE}.log") 2>&1

# 1) Discover unique lanes for this SAMPLE
mapfile -t LANES < <(
  awk -F'\t' -v s="$SAMPLE" '$1==s {
    if ( match($3,/L[0-9]+/) ) print substr($3,RSTART,RLENGTH)
  }' "$MANIFEST" | sort -u
)

declare -a LANE_BAMS=()

for LANE in "${LANES[@]}"; do
  echo "[$(timestamp)] → Processing ${SAMPLE} lane ${LANE}"

  # Grab R1 and R2 for this sample+lane
  R1_PATH=$(awk -F'\t' -v s="$SAMPLE" -v l="$LANE" '$1==s && $2==1 && $3~l {print $3}' "$MANIFEST")
  R2_PATH=$(awk -F'\t' -v s="$SAMPLE" -v l="$LANE" '$1==s && $2==2 && $3~l {print $3}' "$MANIFEST")

  echo "    R1: $R1_PATH"
  echo "    R2: $R2_PATH"

  # sanity check
  if [[ ! -s "$R1_PATH" || ! -s "$R2_PATH" ]]; then
    echo "[$(timestamp)] ❌ Missing FASTQ for ${SAMPLE}_${LANE}, skipping" >&2
    continue
  fi

  OUT_PREF="${OUTDIR}/${SAMPLE}_${LANE}"
  RG="@RG\tID:${SAMPLE}_${LANE}\tSM:${SAMPLE}\tPL:illumina\tLB:${SAMPLE}_lib\tPU:${LANE}"

  #
  # 2) SeqPrep2 trimming & merge
  #
  start=$(date +%s)
  echo "[$(timestamp)] START SeqPrep2 ${SAMPLE}_${LANE}"
  "$SEQPREP" \
    -f "$R1_PATH" -r "$R2_PATH" \
    -1 "${OUT_PREF}_R1_unmerged.fastq.gz" \
    -2 "${OUT_PREF}_R2_unmerged.fastq.gz" \
    -3 "${OUT_PREF}_R1_unpaired.fastq.gz" \
    -4 "${OUT_PREF}_R2_unpaired.fastq.gz" \
    -s "${OUT_PREF}_merged.fastq.gz" \
    -q 15 -L 30 -o 15 -m 0.05 \
    -A AGATCGGAAGAGCACACGTC -B AGATCGGAAGAGCGTCGTGT \
    -C ATCTCGTATGCCGTCTTCTGCTTG \
    -D GATCTCGGTGGTCGCCGTATCATT \
    -S "${OUT_PREF}_seqprep.log"
  zcat "${OUT_PREF}_R1_unpaired.fastq.gz" "${OUT_PREF}_R2_unpaired.fastq.gz" \
    | gzip -c > "${OUT_PREF}_singleton.fastq.gz"
  end=$(date +%s)
  echo "[$(timestamp)] DONE SeqPrep2 (Duration: $((end-start))s)"

  #
  # 3) BWA aln/samse & sort (merged)
  #
  start=$(date +%s)
  echo "[$(timestamp)] START BWA samse (merged)"
  bwa samse "$REFERENCE" \
    <(bwa aln -l 1024 -n 0.03 -o 2 -t "$THREADS" "$REFERENCE" "${OUT_PREF}_merged.fastq.gz") \
    "${OUT_PREF}_merged.fastq.gz" -r "$RG" \
    | samtools sort -@ "$THREADS" -o "${OUT_PREF}_merged_tmp.bam"
  end=$(date +%s)
  echo "[$(timestamp)] DONE samse (merged) (Duration: $((end-start))s)"

  #
  # 4) BWA aln/sampe & sort (unmerged)
  #
  start=$(date +%s)
  echo "[$(timestamp)] START BWA sampe (unmerged)"
  bwa sampe "$REFERENCE" \
    <(bwa aln -l 1024 -n 0.03 -o 2 -t "$THREADS" "$REFERENCE" "${OUT_PREF}_R1_unmerged.fastq.gz") \
    <(bwa aln -l 1024 -n 0.03 -o 2 -t "$THREADS" "$REFERENCE" "${OUT_PREF}_R2_unmerged.fastq.gz") \
    "${OUT_PREF}_R1_unmerged.fastq.gz" "${OUT_PREF}_R2_unmerged.fastq.gz" -r "$RG" \
    | samtools sort -@ "$THREADS" -o "${OUT_PREF}_unmerged_tmp.bam"
  end=$(date +%s)
  echo "[$(timestamp)] DONE sampe (unmerged) (Duration: $((end-start))s)"

  #
  # 5) BWA aln/samse & sort (singleton)
  #
  start=$(date +%s)
  echo "[$(timestamp)] START BWA samse (singleton)"
  bwa samse "$REFERENCE" \
    <(bwa aln -l 1024 -n 0.03 -o 2 -t "$THREADS" "$REFERENCE" "${OUT_PREF}_singleton.fastq.gz") \
    "${OUT_PREF}_singleton.fastq.gz" -r "$RG" \
    | samtools sort -@ "$THREADS" -o "${OUT_PREF}_singleton_tmp.bam"
  end=$(date +%s)
  echo "[$(timestamp)] DONE samse (singleton) (Duration: $((end-start))s)"

  #
  # 6) QC: flagstat
  #
  echo "[$(timestamp)] START flagstat"
  samtools flagstat "${OUT_PREF}_merged_tmp.bam"    > "${OUT_PREF}_flagstat_merged.txt"
  samtools flagstat "${OUT_PREF}_unmerged_tmp.bam"  > "${OUT_PREF}_flagstat_unmerged.txt"
  samtools flagstat "${OUT_PREF}_singleton_tmp.bam" > "${OUT_PREF}_flagstat_singleton.txt"
  echo "[$(timestamp)] DONE flagstat"

  #
  # 7) Merge per‑lane BAMs
  #
  echo "[$(timestamp)] START lane‑merge"
  merge_inputs=()
  for t in merged unmerged singleton; do
    b="${OUT_PREF}_${t}_tmp.bam"
    if [[ -s "$b" && $(samtools view -c "$b") -gt 0 ]]; then
      merge_inputs+=( "$b" )
    else
      echo "[$(timestamp)] ⚠️ Skipping empty $b"
    fi
  done

  if (( ${#merge_inputs[@]} == 0 )); then
    echo "[$(timestamp)] ❌ No reads for ${SAMPLE}_${LANE}" >&2
    exit 1
  fi

  start=$(date +%s)
  samtools merge -@ "$THREADS" - "${merge_inputs[@]}" \
    | samtools sort -@ "$THREADS" -o "${OUTDIR}/${SAMPLE}_${LANE}.lane.bam"
  end=$(date +%s)
  echo "[$(timestamp)] DONE lane‑merge (Duration: $((end-start))s)"

  samtools index "${OUTDIR}/${SAMPLE}_${LANE}.lane.bam"
  rm "${OUT_PREF}"_*_tmp.bam
  LANE_BAMS+=( "${OUTDIR}/${SAMPLE}_${LANE}.lane.bam" )

done

#
# 8) Merge across all lanes
#
echo "[$(timestamp)] START sample‑merge across ${#LANE_BAMS[@]} lane(s)"
start=$(date +%s)
if (( ${#LANE_BAMS[@]} == 1 )); then
  mv "${LANE_BAMS[0]}" "${OUTDIR}/${SAMPLE}.bam"
else
  samtools merge -@ "$THREADS" - "${LANE_BAMS[@]}" \
    | samtools sort -@ "$THREADS" -o "${OUTDIR}/${SAMPLE}.bam"
  rm -f "${LANE_BAMS[@]}"
fi
end=$(date +%s)
echo "[$(timestamp)] DONE sample‑merge (Duration: $((end-start))s)"

samtools index "${OUTDIR}/${SAMPLE}.bam"
samtools flagstat "${OUTDIR}/${SAMPLE}.bam"   > "${OUTDIR}/${SAMPLE}_flagstat.txt"
samtools coverage "${OUTDIR}/${SAMPLE}.bam"   > "${OUTDIR}/${SAMPLE}_coverage.txt"

#
# 9) Mark duplicates
#
echo "[$(timestamp)] START markdup"
start=$(date +%s)
samtools sort -n -@ "$THREADS" \
    -o "${OUTDIR}/${SAMPLE}.namesort.bam" "${OUTDIR}/${SAMPLE}.bam"
samtools fixmate -m -@ "$THREADS" \
    "${OUTDIR}/${SAMPLE}.namesort.bam" "${OUTDIR}/${SAMPLE}.fixmate.bam"
samtools sort -@ "$THREADS" \
    -o "${OUTDIR}/${SAMPLE}.fixmate.sorted.bam" "${OUTDIR}/${SAMPLE}.fixmate.bam"
samtools markdup -r -@ "$THREADS" \
    "${OUTDIR}/${SAMPLE}.fixmate.sorted.bam" "${OUTDIR}/${SAMPLE}.rmdup.bam"
end=$(date +%s)
echo "[$(timestamp)] DONE markdup (Duration: $((end-start))s)"

rm -f "${OUTDIR}/${SAMPLE}".namesort.bam \
      "${OUTDIR}/${SAMPLE}".fixmate.bam \
      "${OUTDIR}/${SAMPLE}".fixmate.sorted.bam

samtools index "${OUTDIR}/${SAMPLE}.rmdup.bam"
samtools flagstat "${OUTDIR}/${SAMPLE}.rmdup.bam" \
    > "${OUTDIR}/${SAMPLE}_rmdup_flagstat.txt"
samtools coverage "${OUTDIR}/${SAMPLE}.rmdup.bam" \
    > "${OUTDIR}/${SAMPLE}_rmdup_coverage.txt"

#
# 10) mapDamage
#
echo "[$(timestamp)] START mapDamage"
start=$(date +%s)
mapDamage -i "${OUTDIR}/${SAMPLE}.rmdup.bam" -r "$REFERENCE" \
    -l 150 \
    -d "${OUTDIR}/out_mapdamage/mapdamage_${SAMPLE}" \
    -y "$MAX_MISINCORP_FREQUENCY" \
    -n "$MAPDOWNSAMPLE" \
    --merge-reference-sequences \
    -m "$READ_PLOT_LENGTH" \
    -t "${REF_NAME}_${SAMPLE}"
end=$(date +%s)
echo "[$(timestamp)] DONE mapDamage (Duration: $((end-start))s)"

echo "[$(timestamp)] ✅ Finished processing ${SAMPLE}"