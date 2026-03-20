#!/bin/bash -l
#SBATCH --job-name=ab_qc
#SBATCH --mail-user=katie.herron@ucdconnect.ie
#SBATCH --mail-type=END
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --time=240:00:00
#SBATCH --array=0-42%8
#SBATCH --output=logs/%x-%A_%a.slurm

set -euo pipefail; mkdir -p logs

# Input dir with per-sample vcfs
VCF_DIR="${1:-/home/people/22211215/scratch/tridubire/03-variants/bcftools_subgen_2n}"
OUT_DIR="$VCF_DIR/ab_qc"
ROOT="/home/people/22211215/scratch"
WD="$ROOT/tridubire"
IMG_BCFTOOLS="docker://quay.io/biocontainers/bcftools:1.22--h3a4d415_1"

mkdir -p "$OUT_DIR"
shopt -s nullglob
mapfile -t VCFS < <(ls -1 "$VCF_DIR"/*.vcf.gz 2>/dev/null | sort)
if (( ${#VCFS[@]} == 0 )); then
  echo "No VCFs found in $VCF_DIR"; exit 0
fi

if (( SLURM_ARRAY_TASK_ID >= ${#VCFS[@]} )); then
  echo "Array index ${SLURM_ARRAY_TASK_ID} exceeds VCF count ${#VCFS[@]} — nothing to do."; exit 0
fi
VCF="${VCFS[$SLURM_ARRAY_TASK_ID]}"
SAMPLEBASE=$(basename "$VCF" .vcf.gz)
TSV_GZ="$OUT_DIR/${SAMPLEBASE}.AB.tsv.gz"
SUM="$OUT_DIR/${SAMPLEBASE}.AB.summary.txt"

echo "Processing: $VCF"

# Extract AD -> AB table
singularity exec --no-home --bind "$ROOT:$ROOT" --pwd "$WD" "$IMG_BCFTOOLS" \
  bash -lc '
    bcftools query -f "[%CHROM\t%POS\t%SAMPLE\t%AD\n]" "'"$VCF"'" \
    | awk -F"\t" '"'"'{
        if ($4 == "." || $4 == "") next;
        n=split($4,a,",");
        ref=(n>=1 && a[1]~/^[0-9]+$/)?a[1]:0;
        alt=(n>=2 && a[2]~/^[0-9]+$/)?a[2]:0;
        tot=ref+alt;
        if (tot>0) {
          ab=alt/tot;
          printf("%s\t%s\t%s\t%d\t%d\t%.6f\n",$1,$2,$3,ref,alt,ab);
        }
      }'"'"' \
    | bgzip -c > "'"$TSV_GZ"'"
    tabix -f "'"$TSV_GZ"'" 2>/dev/null || true
  '
echo "AB table: $TSV_GZ"

# Quick summary
gzip -cd "$TSV_GZ" \
| awk -F"\t" '
    {ab=$6; c++; sum+=ab; a[c]=ab}
    END{
      if(c==0){print "No AB entries"; exit}
      # simple quantiles
      n=c; asort(a);
      q25=a[int(0.25*(n+1))]; q50=a[int(0.50*(n+1))]; q75=a[int(0.75*(n+1))];
      printf("N_ab=%d\nmean=%.4f\np25=%.4f\nmedian=%.4f\np75=%.4f\n",
             c, sum/c, q25, q50, q75);
    }' > "$SUM"
echo "Summary: $SUM"

# Quick hist
if command -v Rscript >/dev/null 2>&1; then
  PDF="${TSV_GZ%.tsv.gz}.AB_hist.pdf"
  Rscript - <<'RS' "$TSV_GZ" "$PDF"
args <- commandArgs(trailingOnly=TRUE)
tsv <- args[1]; pdfout <- args[2]
con <- gzfile(tsv, "rt")
d <- tryCatch(read.table(con, header=FALSE, sep="\t",
                         col.names=c("chrom","pos","sample","ADref","ADalt","AB")),
              finally=close(con))
if (nrow(d)>0) {
  pdf(pdfout, width=6, height=4)
  hist(d$AB, breaks=50, xlab="ALT / (REF+ALT)", ylab="Count",
       main=paste0("Allele balance: ", basename(tsv)))
  abline(v=0.5, lty=2)
  dev.off()
}
RS
  echo "Histogram: $PDF"
else
  echo "Skipped"
fi

echo "Done: $SAMPLEBASE"
