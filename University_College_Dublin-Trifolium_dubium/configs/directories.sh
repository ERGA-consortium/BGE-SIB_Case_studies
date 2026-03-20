# configs/directories.sh
# Project paths

ROOT="${ROOT:-/home/people/22211215/scratch}"
PROJECT_NAME="${PROJECT_NAME:-genuspopgen}"

WD="${WD:-${ROOT}/${PROJECT_NAME}}"

# Project directories
DATA_DIR="${DATA_DIR:-${WD}/00-data}"
QC_DIR="${QC_DIR:-${WD}/01-initial_qc}"
ALIGN_DIR="${ALIGN_DIR:-${WD}/02-alignment}"
VARIANT_DIR="${VARIANT_DIR:-${WD}/03-variant_calling}"
FILTER_DIR="${FILTER_DIR:-${WD}/04-filtering}"
POSTVCFQC_DIR="${POSTQC_DIR:-${WD}/05-post_vcf_qc}"
CONFIG_DIR="${CONFIG_DIR:-${WD}/configs}"
ENV_DIR="${ENV_DIR:-${WD}/envs}"
LOG_DIR="${LOG_DIR:-${WD}/logs}"

# DATA_DIR subdirectories
RAW_DIR="${RAW_DIR:-${DATA_DIR}/raw}"
SRA_CACHE_DIR="${SRA_CACHE_DIR:-${DATA_DIR}/sra_cache}"
TMP_DIR="${TMP_DIR:-${DATA_DIR}/tmp}"
FQDUMP_TMP="${FQDUMP_TMP:-${TMP_DIR}/fasterq}"

# QC_DIR subdirectories
FQS_CONF_DIR="${FQS_CONF_DIR:-${QC_DIR}/fastq_screen}"
FQS_DB_DIR="${FQS_DB_DIR:-${FQS_CONF_DIR}/db}"
FQS_OUT_DIR="${FQS_OUT_DIR:-${QC_DIR}/fastq_screen_out}"

RAW_FASTQC_DIR="${RAW_FASTQC_DIR:-${QC_DIR}/fastqc"
TRIM_DIR="${TRIM_DIR:-${QC_DIR}/fastp}"

# ALIGN_DIR subdirectories
RUN_BAM_DIR="${RUN_BAM_DIR:-${ALIGN_DIR}/bam-per_run}"
SAMPLE_BAM_DIR="${SAMPLE_BAM_DIR:-${ALIGN_DIR}/bam-merged}"
FLAGSTAT_DIR="${FLAGSTAT_DIR:-${ALIGN_DIR}/flagstat}"
IDXSTAT_DIR="${IDXSTAT_DIR:-${ALIGN_DIR}/idxstat}"

# Create directories
mkdir -p "$DATA_DIR" "$QC_DIR" "$ALIGN_DIR" "$VARIANT_DIR" "$FILTER_DIR" "$POSTVCFQC_DIR" "$FQS_CONF_DIR" "$FQS_DB_DIR" "$FQS_OUT_DIR" "$RAW_DIR" "$SRA_CACHE_DIR" "$TMP_DIR" "$FQDUMP_TMP" "$RAW_FASTQC_DIR" "$TRIM_DIR" "$CONFIG_DIR" "$ENV_DIR" "$LOG_DIR" "$RUN_BAM_DIR" "$SAMPLE_BAM_DIR" "$FLAGSTAT_DIR" "$IDXSTAT_DIR"

# Project-managed files
SAMPLESHEET="${SAMPLESHEET:-${CONFIG_DIR}/samplesheet.csv}"
FULL_REF_FASTA="${FULL_REF_FASTA:-${DATA_DIR}/triDubi.fna}"
REF_FASTA="${REF_FASTA:-${DATA_DIR}/triDubi_chromosomes.fna}"
FQS_CONF="${FQS_CONF:-${FQS_CONF_DIR}/fastq_screen.conf}"
