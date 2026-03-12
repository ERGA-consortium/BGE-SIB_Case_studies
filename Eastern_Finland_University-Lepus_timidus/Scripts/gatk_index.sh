#!/bin/bash

#SBATCH --job-name=hyr_gatkindex
#SBATCH --account=project_2002674
#SBATCH --time=15:00
#SBATCH --mem-per-cpu=2G
#SBATCH --cpus-per-task=4
#SBATCH --partition=test
#SBATCH --output=HyR_gatkindex_out_%j.txt
#SBATCH --error=HyR_gatkindex_err_%j.txt
#SBATCH --mail-type=END

#variables
refdir=/scratch/project_2002674/zsofia/Hybrid_rerun/GENOMES/

#modules
export PATH="/projappl/project_2002674/chromcomp_env//bin:$PATH"

#Run
for i in $refdir/*.fa; do
    gatk CreateSequenceDictionary --java-options "-Xmx8G" -R $i;
done
