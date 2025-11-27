#!/bin/bash -l

#SBATCH -A hpc2n2025-120	
#SBATCH -n 1
#SBATCH -c 8
#SBATCH -t 1-00:00:00
#SBATCH -J qc
#SBATCH --output=/proj/nobackup/hpc2nstor2025-059/mimmi/mind_the_gaps/reports/sbatch/qc/sbatch_R-%x_%j.out
#SBATCH --error=/proj/nobackup/hpc2nstor2025-059/mimmi/mind_the_gaps/reports/sbatch/qc/sbatch_R-%x_%j.err

species=$1
echo ${species}

eval "$(mamba shell hook --shell bash)"
mamba activate chromcomp

cd /proj/nobackup/hpc2nstor2025-059/mimmi/mind_the_gaps

fastqc -t 8 data/${species}/reads/raw/*.fq.gz -o data/${species}/qc/fastqc

multiqc data/${species}/qc/fastqc -o data/${species}/qc