#!/bin/bash -l

#SBATCH -A hpc2n2025-120	
#SBATCH -n 1
#SBATCH -c 1
#SBATCH -t 1-00:00:00
#SBATCH -J indexing
#SBATCH --output=/proj/nobackup/hpc2nstor2025-059/mimmi/mind_the_gaps/reports/sbatch/indexing/sbatch_R-%x_%j-%a.out
#SBATCH --error=/proj/nobackup/hpc2nstor2025-059/mimmi/mind_the_gaps/reports/sbatch/indexing/sbatch_R-%x_%j-%a.err
#SBATCH -a 0-1

eval "$(mamba shell hook --shell bash)"
mamba activate chromcomp

cd /proj/nobackup/hpc2nstor2025-059/mimmi/mind_the_gaps

species=$1

run_list=(draft chrom)
asm=$(echo ${run_list[$SLURM_ARRAY_TASK_ID]})

echo ${species}
echo ${asm}

# Build FM‑index for BWA MEM (necessary once per assembly).
bwa index data/${species}/assemblies/${asm}/${asm}_round-2_masked.fa

# Generate FASTA index (.fai) so downstream tools know chromosome lengths.
samtools faidx data/${species}/assemblies/${asm}/${asm}_round-2_masked.fa

# Create a sequence dictionary (.dict) required by GATK.
gatk CreateSequenceDictionary -R data/${species}/assemblies/${asm}/${asm}_round-2_masked.fa
