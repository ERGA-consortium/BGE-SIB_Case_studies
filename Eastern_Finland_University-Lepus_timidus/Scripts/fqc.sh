#!/bin/bash

#SBATCH --job-name=MtG_fqc
#SBATCH --account=project_2002674
#SBATCH --time=6:00:00
#SBATCH --mem-per-cpu=100M
#SBATCH --cpus-per-task=8
#SBATCH --partition=small
#SBATCH --mail-type=END

#SBATCH --output=MtG_fqc_output_%j.txt
#SBATCH --error=MtG_fqc_errors_%j.txt


#Process_radtags
module load biokit
export PATH="/projappl/project_2002674/fqc-0.12.1/FastQC/:$PATH"

#Variables
fqdir=/scratch/project_2002674/zsofia/Mind_the_Gap/00_Data/clean_demulti/
outdir=/scratch/project_2002674/zsofia/Mind_the_Gap/01_QC/demultiplexed/

srun fastqc -o $outdir/LE/ $fqdir/LE/*.fq.gz
srun fastqc -o $outdir/LT/ $fqdir/LT/*.fq.gz
srun fastqc -o $outdir/LT_parental/ $fqdir/LT_parental/*.fq.gz
srun fastqc -o $outdir/LE_parental/ $fqdir/LE_parental/*.fq.gz
