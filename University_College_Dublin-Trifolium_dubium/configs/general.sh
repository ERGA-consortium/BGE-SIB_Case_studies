# configs/general.sh
# Run settings

PLOIDY="${PLOIDY:-2}"
THREADS="${THREADS:-8}"
MEM_GB="${MEM_GB:-32}"

# Singularity run settings - may differ depending on HPC/local set-up
SINGULARITY_EXEC="${SINGULARITY_EXEC:-singularity exec}"
SING_EXTRA_ARGS="${SING_EXTRA_ARGS:---no-home}"
SING_BIND="${SING_BIND:-${ROOT}:${ROOT}}"
SING_PWD="${SING_PWD:-${WD}}"

CONTAINER_MODE="${CONTAINER_MODE:-enabled}"  # or "disabled" if not using containers, have all software installed/present in a conda env
