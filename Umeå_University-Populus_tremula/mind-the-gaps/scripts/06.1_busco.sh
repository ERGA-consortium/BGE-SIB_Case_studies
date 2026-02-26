#!/bin/bash -l

#SBATCH -A hpc2n2025-120	
#SBATCH -n 1
#SBATCH -c 80
##SBATCH -C zen4
#SBATCH -t 2-00:00:00
#SBATCH -J busco
#SBATCH --output=/proj/nobackup/hpc2nstor2025-059/mimmi/mind_the_gaps/reports/sbatch/busco/sbatch_R-%x_%j.out
#SBATCH --error=/proj/nobackup/hpc2nstor2025-059/mimmi/mind_the_gaps/reports/sbatch/busco/sbatch_R-%x_%j.err


eval "$(mamba shell hook --shell bash)"
mamba activate chromcomp

cd /proj/nobackup/hpc2nstor2025-059/mimmi/mind_the_gaps

# run_list=(draft chrom)
# asm=$(echo ${run_list[$SLURM_ARRAY_TASK_ID]})
asm="draft"
species=$1

echo ${species}
echo ${asm}

busco -i  data/${species}/assemblies/${asm}/${asm}_round-2_masked.fa \
  --mode genome \
  -f \
  --cpu 2 \
  --lineage_dataset /proj/nobackup/hpc2nstor2025-059/mimmi/mind_the_gaps/data/busco_db/viridiplantae_odb12 \
  -o data/${species}/annotation/busco/${asm}/${asm}_busco
