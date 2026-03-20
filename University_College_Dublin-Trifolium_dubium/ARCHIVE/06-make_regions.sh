#!/bin/bash -l
#SBATCH --job-name=regions
#SBATCH --mail-user=katie.herron@ucdconnect.ie
#SBATCH --mail-type=END
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --time=12:00:00
#SBATCH --output=logs/%x-%A_%a.slurm

# Paths
ROOT="/home/people/22211215/scratch"
WD="$ROOT/tridubire"
REF="$WD/00-data/triDubi-GCA_951804385.1.fna"
REGIONS_DIR="$WD/00-data/regions"

# Container
IMG_SAMTOOLS="docker://quay.io/biocontainers/samtools:1.22.1--h96c455f_0"

mkdir -p "$REGIONS_DIR"

if [[ ! -s "${REF}.fai" ]]; then
  echo "Indexing reference..."
  singularity exec --no-home --bind "$ROOT:$ROOT" --pwd "$WD" "$IMG_SAMTOOLS" \
    samtools faidx "$REF"
else
  echo "Found existing ${REF}.fai"
fi

awk 'BEGIN{OFS="\t"}{print $1,0,$2}' "${REF}.fai" > "$REGIONS_DIR/triDubi-genome.bed"

cat > "$REGIONS_DIR/triDubi-SG1.contigs.txt" <<'EOF'
OX638062.1
OX638063.1
OX638064.1
OX638065.1
OX638068.1
OX638069.1
OX638071.1
EOF

cat > "$REGIONS_DIR/triDubi-SG2.contigs.txt" <<'EOF'
OX638066.1
OX638067.1
OX638070.1
OX638072.1
OX638073.1
OX638074.1
OX638075.1
OX638076.1
EOF

cat > "$REGIONS_DIR/triDubi-mito.contigs.txt" <<'EOF'
OX638077.1
OX638078.1
EOF

cat > "$REGIONS_DIR/triDubi-chloro.contigs.txt" <<'EOF'
OX638079.1
EOF

cat "$REGIONS_DIR/triDubi-mito.contigs.txt" "$REGIONS_DIR/triDubi-chloro.contigs.txt" \
  | sort -u > "$REGIONS_DIR/triDubi-org.contigs.txt"

cat "$REGIONS_DIR/triDubi-SG1.contigs.txt" "$REGIONS_DIR/triDubi-SG2.contigs.txt" \
  | sort -u > "$REGIONS_DIR/triDubi-nuclear.contigs.txt"

# Unplaced scaffolds
cat "$REGIONS_DIR/triDubi-SG1.contigs.txt" "$REGIONS_DIR/triDubi-SG2.contigs.txt" "$REGIONS_DIR/triDubi-org.contigs.txt" \
  | sort -u > "$REGIONS_DIR/triDubi-used.contigs.txt"
cut -f1 "${REF}.fai" | grep -vxF -f "$REGIONS_DIR/triDubi-used.contigs.txt" > "$REGIONS_DIR/triDubi-unplaced.contigs.txt" || true

# Builds beds
mkbed() {
  local list="$1" bed="$2"
  if [[ -s "$list" ]]; then
    grep -Ff "$list" "$REGIONS_DIR/triDubi-genome.bed" > "$bed" || true
  else
    : > "$bed"
  fi
}

mkbed "$REGIONS_DIR/triDubi-SG1.contigs.txt"      "$REGIONS_DIR/triDubi-SG1.bed"
mkbed "$REGIONS_DIR/triDubi-SG2.contigs.txt"      "$REGIONS_DIR/triDubi-SG2.bed"
mkbed "$REGIONS_DIR/triDubi-nuclear.contigs.txt"  "$REGIONS_DIR/triDubi-nuclear.bed"
mkbed "$REGIONS_DIR/triDubi-mito.contigs.txt"     "$REGIONS_DIR/triDubi-mito.bed"
mkbed "$REGIONS_DIR/triDubi-chloro.contigs.txt"   "$REGIONS_DIR/triDubi-chloro.bed"
mkbed "$REGIONS_DIR/triDubi-org.contigs.txt"      "$REGIONS_DIR/triDubi-org.bed"
mkbed "$REGIONS_DIR/triDubi-unplaced.contigs.txt" "$REGIONS_DIR/triDubi-unplaced.bed"

echo "BEDs written to: $REGIONS_DIR"
ls -1 "$REGIONS_DIR"/triDubi-*.bed || true
