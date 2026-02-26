#!/bin/bash -l

#SBATCH -A hpc2n2025-120	
#SBATCH -n 1
#SBATCH -c 1
#SBATCH -t 1-00:00:00
#SBATCH -J trimming
#SBATCH --output=/proj/nobackup/hpc2nstor2025-059/mimmi/mind_the_gaps/reports/sbatch/trimming/sbatch_R-%x_%j-%a.out
#SBATCH --error=/proj/nobackup/hpc2nstor2025-059/mimmi/mind_the_gaps/reports/sbatch/trimming/sbatch_R-%x_%j-%a.err
#SBATCH -a 0-9

eval "$(mamba shell hook --shell bash)"
mamba activate chromcomp

cd /proj/nobackup/hpc2nstor2025-059/mimmi/mind_the_gaps

species=$1

sample_list=($(cut -f1 data/${species}/sample_info.txt))
sample_id=$(echo ${sample_list[${SLURM_ARRAY_TASK_ID}]})

echo ${species}
echo ${sample_id}

### quality filter with fastp
# -i forward reads
# -I reverse reads
# -o trimmed output (R1)
# -O trimmed output (R2)
# -q minimum mean quality per sliding window
# -u discard read if >30 % bases fall below Q20
# -l minimum post‑trim length (bp)
# -g polyG tail trimming for Illumina NextSeq/NovaSeq data
# -h  interactive report
# -j  machine‑parsable stats
fastp \
    -i data/${species}/reads/raw/${sample_id}_R1.fq.gz \
    -I data/${species}/reads/raw/${sample_id}_R2.fq.gz \
    -o data/${species}/reads/trimmed/${sample_id}_R1.trim.fq.gz \
    -O data/${species}/reads/trimmed/${sample_id}_R2.trim.fq.gz \
    -q 20 \
    -u 30 \
    -l 50 \
    -g \
    -h data/${species}/qc/fastp/html/${sample_id}.html \
    -j data/${species}/qc/fastp/json/${sample_id}.json