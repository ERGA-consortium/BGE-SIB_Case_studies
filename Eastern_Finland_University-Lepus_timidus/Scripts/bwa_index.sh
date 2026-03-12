#!/bin/bash

#SBATCH --job-name=hyr_bwaindex
#SBATCH --account=project_2002674
#SBATCH --time=2:00:00
#SBATCH --mem-per-cpu=5G
#SBATCH --cpus-per-task=5
#SBATCH --partition=small
#SBATCH --output=HyR_index_out_%j.txt
#SBATCH --error=HyR_index_err_%j.txt
#SBATCH --mail-type=END

#variables
refdir=/scratch/project_2002674/zsofia/Hybrid_hares/REFERENCE/

#modules
module load bwa-mem2/2.2.1

#Run
for i in $refdir/*.fa; do
    bwa index $i;
done
