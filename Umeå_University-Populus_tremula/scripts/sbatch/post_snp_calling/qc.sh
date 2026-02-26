#!/bin/bash -l

#SBATCH -A hpc2n	
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 8
#SBATCH -t 1-00:00:00
#SBATCH -J qc
#SBATCH --output=reports/sbatch/qc/sbatch_R-%x_%j.out
#SBATCH --error=reports/sbatch/qc/sbatch_R-%x_%j.err

ml GCC/13.2.0 
ml OpenMPI/4.1.6
ml FastQC/0.12.1-Java-11
ml MultiQC/1.22.3

cd 

fq_dir="fastq/NordAsp/raw/reads/"

fastqc -t 8 ${fq_dir}*.fq.gz -o output_data/qc/fastqc

multiqc output_data/qc/fastqc -o output_data/qc