#!/bin/bash -l
#SBATCH --job-name=summarise_alignment_qc
#SBATCH --cpus-per-task=2
#SBATCH --time=02:00:00
#SBATCH --output=logs/%x-%A.slurm

ROOT="/home/people/22211215/scratch/"
WD="$ROOT/tridubire"
ALIGN_DIR="$WD/02-alignment"
FLAGSTAT_DIR="$ALIGN_DIR/flagstat"
IDXSTAT_DIR="$ALIGN_DIR/idxstat"
MOSDEPTH_DIR="$ALIGN_DIR/mosdepth"
OUT_DIR="$ALIGN_DIR/qc_summary"
mkdir -p "$OUT_DIR"

OUT_TSV="$OUT_DIR/qc_summary.tsv"
printf "sample\ttotal_reads\tmapped_reads\tmapped_rate\tproperly_paired_reads\tproperly_paired_rate\tduplicate_reads\tmean_depth\tcovered_fraction\n" > "$OUT_TSV"

shopt -s nullglob
for FS in "$FLAGSTAT_DIR"/*.markdup.flagstat.txt; do
  sample="$(basename "$FS" .markdup.flagstat.txt)"

  MS="$MOSDEPTH_DIR/${sample}.mosdepth.summary.txt"
  GD="$MOSDEPTH_DIR/${sample}.mosdepth.global.dist.txt"

  # --- Parse flagstat robustly ---
  total_reads=$(awk 'match($0,/^([0-9]+) \+ [0-9]+ in total/,a){print a[1]; exit}' "$FS"); [[ -z "$total_reads" ]] && total_reads="NA"
  mapped_reads=$(awk 'match($0,/^([0-9]+) \+ [0-9]+ mapped/,a){print a[1]; exit}' "$FS"); [[ -z "$mapped_reads" ]] && mapped_reads="NA"
  mapped_rate=$(grep -m1 ' mapped ' "$FS" | awk 'match($0,/\(([0-9.]+)%/,a){print a[1]}' ); [[ -z "$mapped_rate" ]] && mapped_rate="NA"
  properly_paired_reads=$(awk 'match($0,/^([0-9]+) \+ [0-9]+ properly paired/,a){print a[1]; exit}' "$FS"); [[ -z "$properly_paired_reads" ]] && properly_paired_reads="NA"
  properly_paired_rate=$(grep -m1 ' properly paired ' "$FS" | awk 'match($0,/\(([0-9.]+)%/,a){print a[1]}' ); [[ -z "$properly_paired_rate" ]] && properly_paired_rate="NA"
  duplicate_reads=$(awk 'match($0,/^([0-9]+) \+ [0-9]+ duplicates$/,a){print a[1]; exit}' "$FS"); [[ -z "$duplicate_reads" ]] && duplicate_reads="NA"

  # Mosdepth mean depth
  if [[ -s "$MS" ]]; then
    mean_depth=$(awk '$1=="total"{print $4}' "$MS")
  else
    mean_depth="NA"
  fi
  
  if [[ -s "$GD" ]]; then
    cov1=$(awk '$1=="total" && $2==1 {print $3; found=1} END{if(!found) print ""}' "$GD")
    if [[ -n "$cov1" ]]; then
      covered_fraction="$cov1"
    else
      frac0=$(awk '$1=="total" && $2==0 {print $3; found=1} END{if(!found) print ""}' "$GD")
      if [[ -n "$frac0" ]]; then
        covered_fraction=$(awk -v f="$frac0" 'BEGIN{printf "%.6f", 1.0 - f}')
      else
        covered_fraction="NA"
      fi
    fi
    
  # fallback to region.dist
  elif [[ -s "$RD" && -s "$MS" ]]; then
    covered_fraction=$(awk '
      FNR==NR && $1!="total" {len[$1]=$2; tot+=$2; next}
      FNR!=NR && $2==1 && ($1 in len) {cov += len[$1]*$3}
      END{ if(tot>0) printf "%.6f\n", cov/tot; else print "NA" }
    ' "$MS" "$RD")
  else
    covered_fraction="NA"
  fi

  printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
    "$sample" "$total_reads" "$mapped_reads" "$mapped_rate" "$properly_paired_reads" "$properly_paired_rate" "$duplicate_reads" "$mean_depth" "$covered_fraction" \
    >> "$OUT_TSV"
done

echo "Wrote: $OUT_TSV"
