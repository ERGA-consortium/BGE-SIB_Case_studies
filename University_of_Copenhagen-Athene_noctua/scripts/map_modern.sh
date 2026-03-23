#!/bin/bash
#SBATCH --job-name=modern_bwa
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.err
#SBATCH --cpus-per-task=8
#SBATCH --mem=14G
#SBATCH --time=00:20:00
set -euo pipefail
trap 'echo "❌ Error occurred during processing of ${ID}. Check the log file at ${OUTdir}/pipeline_${ID}.log" >&2' ERR

# Modules
module load mapdamage2/2.2.2
module load gcc R/3.6.1
module load fastp/0.24.0
module load perl openjdk fastqc/0.12.1

# Check args
if [[ $# -ne 5 ]]; then
    echo "Usage: $0 <reference.fasta> <manifest.tsv> <sample_id> <threads> <output_dir>"
    exit 1
fi

REFERENCE=$1
MANIFEST=$2
ID=$3
THREADS=$4
OUTdir=$5

# Derived
mkdir -p "$OUTdir" tmp trimmed out_fastqc out_fastqc_trimmed
PICARD=/opt/software/picard/2.27.5/picard.jar
REF_NAME=$(basename "$REFERENCE" .fasta)
mapDamage_downSample=10000
READ_PLOT_LENGTH=50
MAX_MISINCORP_FREQUENCY=0.4

exec > >(tee -i "${OUTdir}/pipeline_${ID}.log")
exec 2>&1

# Get the unique lanes for ID
lanes=($(awk -v id="$ID" '$1==id && $2==1 {print $4}' "$MANIFEST" | sort | uniq))

for LANE in "${lanes[@]}"; do
    echo "🔧 Processing sample $ID, lane $LANE $(date)"

R1=$(awk -v id="$ID" -v lane="$LANE" '$1==id && $2==1 && $4==lane {print $3; exit}' "$MANIFEST")
R2=$(awk -v id="$ID" -v lane="$LANE" '$1==id && $2==2 && $4==lane {print $3; exit}' "$MANIFEST")

    if [[ -z "$R1" || -z "$R2" ]]; then
        echo "❗ Missing R1 or R2 for $ID lane $LANE"
        exit 1
    fi

    OUT="${ID}_${LANE}"
    TRIM_R1="trimmed/${OUT}_R1.fastq.gz"
    TRIM_R2="trimmed/${OUT}_R2.fastq.gz"

    echo "🚀 Running fastp for $OUT $(date)"
    fastp -w "$THREADS" \
        -i "$R1" -I "$R2" \
        -o "$TRIM_R1" -O "$TRIM_R2" \
        -h "trimmed/${OUT}.html" -j "trimmed/${OUT}.json"

    fastqc "$TRIM_R1" --outdir=out_fastqc_trimmed --threads 4
    fastqc "$TRIM_R2" --outdir=out_fastqc_trimmed --threads 4

    RG="@RG\tID:${OUT}\tSM:${ID}\tPL:illumina\tLB:${ID}_lib\tPU:${LANE}"

	# Choose how many threads samtools sort should get (half, but at least 1)
	SORT_THREADS=$(( THREADS / 2 ))
	(( SORT_THREADS < 1 )) && SORT_THREADS=1

	# Give remaining threads to bwa (at least 1)
	BWA_THREADS=$(( THREADS - SORT_THREADS ))
	(( BWA_THREADS < 1 )) && BWA_THREADS=1

	# samtools sort memory *per thread* (tune this)
	# If you request 14G total and SORT_THREADS=4, then 2G/thread is ~8G total to sort.
	SORT_MEM_PER_THREAD="2G"

	# Temp dir for sort (ideally node-local scratch if available)
	# Common: $SLURM_TMPDIR. Fallback to outdir tmp.
	SORT_TMPDIR="${SLURM_TMPDIR:-${OUTdir}/tmp}"
	mkdir -p "$SORT_TMPDIR"

	echo "🧬 Mapping $OUT $(date)"
	echo "    bwa threads: $BWA_THREADS | sort threads: $SORT_THREADS | sort mem/thread: $SORT_MEM_PER_THREAD | sort tmp: $SORT_TMPDIR"

	bwa mem -t "$BWA_THREADS" -R "$RG" "$REFERENCE" "$TRIM_R1" "$TRIM_R2" \
	  | samtools sort \
		  -@ "$SORT_THREADS" \
		  -m "$SORT_MEM_PER_THREAD" \
		  -T "${SORT_TMPDIR}/${OUT}.sorttmp" \
		  -o "${OUTdir}/${OUT}_tmp.bam" \
		  -

    samtools index "${OUTdir}/${OUT}_tmp.bam"
    samtools flagstat "${OUTdir}/${OUT}_tmp.bam" > "${OUTdir}/flagstat_${OUT}_tmp.bam.txt"
    samtools coverage "${OUTdir}/${OUT}_tmp.bam" > "${OUTdir}/coverage_${OUT}_tmp.bam.txt"
done

echo "📦 Merging all lanes for $ID $(date)"
shopt -s nullglob
tmp_bams=( "${OUTdir}/${ID}"_*_tmp.bam )
shopt -u nullglob

if [[ ${#tmp_bams[@]} -eq 0 ]]; then
  echo "❌ ERROR: No tmp BAMs found for $ID in $OUTdir"
  exit 1
elif [[ ${#tmp_bams[@]} -eq 1 ]]; then
  echo "➡️ Only one tmp BAM for $ID — skipping merge $(date)"
  mv -f "${tmp_bams[0]}" "${OUTdir}/${ID}.bam"
else
  echo "🔗 Merging ${#tmp_bams[@]} BAM files for $ID $(date)"
  # merge to stdout, then sort to final bam
  samtools merge -u -@ "$THREADS" - "${tmp_bams[@]}" \
    | samtools sort -@ "$THREADS" -o "${OUTdir}/${ID}.bam"

  # clean lane BAMs + indexes if merge succeeded
  rm -f "${tmp_bams[@]}" "${tmp_bams[@]/.bam/.bam.bai}"
fi

samtools index "${OUTdir}/${ID}.bam"
samtools flagstat "${OUTdir}/${ID}.bam" > "${OUTdir}/flagstat_${ID}.merged.txt"
samtools coverage "${OUTdir}/${ID}.bam" > "${OUTdir}/coverage_${ID}.merged.txt"

echo "🧹 Marking and removing duplicates $(date)"
samtools sort -n -@ "$THREADS" -o "${OUTdir}/${ID}.namesort.bam" "${OUTdir}/${ID}.bam"
samtools fixmate -m -@ "$THREADS" "${OUTdir}/${ID}.namesort.bam" "${OUTdir}/${ID}.fixmate.bam"
samtools sort -@ "$THREADS" -o "${OUTdir}/${ID}.fixmate.sorted.bam" "${OUTdir}/${ID}.fixmate.bam"

if samtools markdup -r -@ "$THREADS" "${OUTdir}/${ID}.fixmate.sorted.bam" "${OUTdir}/${ID}.rmdup.bam"; then
    echo "✅ Duplicate removal successful for $ID — cleaning up intermediates $(date)"
    rm "${OUTdir}/${ID}.namesort.bam" \
       "${OUTdir}/${ID}.fixmate.bam" \
       "${OUTdir}/${ID}.fixmate.sorted.bam"
else
    echo "❌ ERROR: Duplicate removal failed for $ID. Intermediate files retained for debugging. $(date)"
    exit 1
fi

samtools index "${OUTdir}/${ID}.rmdup.bam"
samtools flagstat "${OUTdir}/${ID}.rmdup.bam" > "${OUTdir}/flagstat_${ID}.rmdup.txt"
samtools coverage "${OUTdir}/${ID}.rmdup.bam" > "${OUTdir}/coverage_${ID}.rmdup.txt"

echo "📉 Running mapDamage $(date)"
mkdir -p out_mapdamage
mapDamage -i "${OUTdir}/${ID}.rmdup.bam" -r "$REFERENCE" \
    -l 150 -d "out_mapdamage/mapdamage_${ID}" \
    -y "$MAX_MISINCORP_FREQUENCY" -n "$mapDamage_downSample" \
    --merge-reference-sequences -m "$READ_PLOT_LENGTH" \
    -t "${REF_NAME}_${ID}"

echo "✅ DONE processing ${ID} $(date)"
echo "--------------------------------------------------"
echo "Final BAM (deduplicated): ${OUTdir}/${ID}.rmdup.bam"
echo "Index:                    ${OUTdir}/${ID}.rmdup.bam.bai"
echo "mapDamage output:         out_mapdamage/mapdamage_${ID}/"
echo  "$(date)"
echo "--------------------------------------------------"