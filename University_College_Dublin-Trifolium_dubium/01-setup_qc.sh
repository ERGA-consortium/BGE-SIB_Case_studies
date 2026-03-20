#!/bin/bash -l
# Project setup

ROOT="${ROOT:-/home/people/22211215/scratch}"
PROJECT_NAME="${PROJECT_NAME:-genuspopgen}"
WD="${WD:-${ROOT}/${PROJECT_NAME}}"

source "$WD/configs/directories.sh"
source "$WD/configs/containers.sh"
source "$WD/configs/general.sh"
source "$WD/bin/lib.sh"

log "Project root:      $WD"
log "Data dir:          $DATA_DIR"
log "Raw reads dir:     $RAW_DIR"
log "Samplesheet:       $SAMPLESHEET"
log "Full ref:          $FULL_REF_FASTA"
log "Chrom-only:        $REF_FASTA"

cd "$FQS_DB_DIR"

# FastQ Screen references
if [[ ! -f human.fna ]]; then
  log "Downloading human (GRCh38) ref for FastQ Screen..."
  wget -O human.fna.gz \
    "ftp://ftp.ensembl.org/pub/release-111/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"
  gunzip -f human.fna.gz
fi

if [[ ! -f ecoli.fna ]]; then
  log "Downloading E. coli (ASM584v2) ref..."
  wget -O ecoli.fna.gz \
    "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.fna.gz"
  gunzip -f ecoli.fna.gz
fi

if [[ ! -f phix174.fna ]]; then
  log "Downloading PhiX174 ref..."
  wget -O phix174.fna.gz \
    "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/819/615/GCF_000819615.1_ViralProj14015/GCF_000819615.1_ViralProj14015_genomic.fna.gz"
  gunzip -f phix174.fna.gz
fi

if [[ ! -f fungi.fna ]]; then
  log "Downloading fungi ref..."
  wget -O fungi.fna.gz \
    "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/182/965/GCF_000182965.3_ASM18296v3/GCF_000182965.3_ASM18296v3_genomic.fna.gz"
  gunzip -f fungi.fna.gz
fi

if [[ ! -f medLupu.fna ]]; then
  log "Downloading Medicago lupulina (drMedLupu1.1) ref as medLupu..."
  wget -O medLupu.fna.gz \
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/958/299/785/GCA_958299785.1_drMedLupu1.1/GCA_958299785.1_drMedLupu1.1_genomic.fna.gz"
  gunzip -f medLupu.fna.gz
fi

if [[ ! -f triDubi_plastid.fna ]]; then
  log "Downloading triDubi plastid genome..."
  wget -O triDubi_plastid.fna \
    "https://www.ncbi.nlm.nih.gov/sviewer/viewer.fcgi?id=OX638079.1&db=nuccore&report=fasta&conwithfeat=on&withparts=on"
fi

if [[ ! -f triDubi_mito1.fna ]]; then
  log "Downloading triDubi mito genome 1..."
  wget -O triDubi_mito1.fna \
    "https://www.ncbi.nlm.nih.gov/sviewer/viewer.fcgi?id=OX638077.1&db=nuccore&report=fasta&conwithfeat=on&withparts=on"
fi

if [[ ! -f triDubi_mito2.fna ]]; then
  log "Downloading triDubi mito genome 2..."
  wget -O triDubi_mito2.fna \
    "https://www.ncbi.nlm.nih.gov/sviewer/viewer.fcgi?id=OX638078.1&db=nuccore&report=fasta&conwithfeat=on&withparts=on"
fi

cd "$WD"

# Download T. dubium reference genome
if [[ ! -f "$FULL_REF_FASTA" ]]; then
  log "Downloading drTriDubi3.1 genome to $FULL_REF_FASTA..."
  tmp="${FULL_REF_FASTA}.gz"
  wget -O "$tmp" \
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/951/804/385/GCA_951804385.1_drTriDubi3.1/GCA_951804385.1_drTriDubi3.1_genomic.fna.gz"
  gunzip -f "$tmp"
else
  log "Skipping: Full triDubi ref already exists at $FULL_REF_FASTA."
fi

if [[ ! -f "$FULL_REF_FASTA" ]]; then
  echo "ERROR: Full triDubi ref not found at $FULL_REF_FASTA" >&2
  exit 1
fi

# Grep headers containing "chrom", for chrom only ref
if [[ ! -f "$REF_FASTA" ]]; then
  log "Extracting chromosome-only triDubi with seqkit"
  run_container "$IMG_SEQKIT" \
    seqkit grep -r -i -n -p "chrom" "$FULL_REF_FASTA" > "$REF_FASTA"
else
  log "Skipping: Chromosome-only triDubi ref already exists."
fi

# BWA indexes
for fasta in "$FQS_DB_DIR"/*.fna; do
  [[ -f "$fasta" ]] || continue

  if [[ -f "${fasta}.bwt.2bit.64" ]]; then
    log "Skipping bwa-mem2 index (exists): $(basename "$fasta")"
    continue
  fi

  log "Indexing FastQ Screen DB with bwa-mem2: $(basename "$fasta")"
  run_container "$IMG_BWAMEM2" bwa-mem2 index "$fasta"
done

for fasta in "$FULL_REF_FASTA" "$REF_FASTA"; do
  [[ -f "$fasta" ]] || continue

  if [[ -f "${fasta}.bwt.2bit.64" ]]; then
    log "Skipping bwa-mem2 index (exists): $(basename "$fasta")"
    continue
  fi

  log "Indexing triDubi ref with bwa-mem2: $(basename "$fasta")"
  run_container "$IMG_BWAMEM2" bwa-mem2 index "$fasta"
done

# samtools faidx indexes
for fasta in "$FULL_REF_FASTA" "$REF_FASTA"; do
  [[ -f "$fasta" ]] || continue

  if [[ -f "${fasta}.fai" ]]; then
    log "Skipping samtools faidx (exists): $(basename "$fasta")"
    continue
  fi

  log "Indexing triDubi ref with samtools faidx: $(basename "$fasta")"
  run_container "$IMG_SAMTOOLS" samtools faidx "$fasta"
done

# Write fastq_screen.conf
log "Writing FastQ Screen config: $FQS_CONF"

{
  echo "# FastQ Screen config"
  echo ""
  echo "DATABASE\thuman_GRCh38\t${FQS_DB_DIR}/human.fna"
  echo "DATABASE\tecoli_K12\t${FQS_DB_DIR}/ecoli.fna"
  echo "DATABASE\tphix174\t${FQS_DB_DIR}/phix174.fna"
  echo "DATABASE\tfungi_ASM18296v3\t${FQS_DB_DIR}/fungi.fna"
  echo "DATABASE\tmedLupu_drMedLupu1.1\t${FQS_DB_DIR}/medLupu.fna"
  echo "DATABASE\ttriDubi_drTriDubi3.1\t${FULL_REF_FASTA}"
  echo "DATABASE\ttriDubi_chromosomes\t${REF_FASTA}"
  echo "DATABASE\ttriDubi_plastid\t${FQS_DB_DIR}/triDubi_plastid.fna"
  echo "DATABASE\ttriDubi_mito1\t${FQS_DB_DIR}/triDubi_mito1.fna"
  echo "DATABASE\ttriDubi_mito2\t${FQS_DB_DIR}/triDubi_mito2.fna"
} > "$FQS_CONF"

log "FastQ Screen config written to: $FQS_CONF"

log "Checking SRA reads listed in samplesheet..."
populate_samplesheet_reads_from_sra "$SAMPLESHEET"

log "Setup complete."
