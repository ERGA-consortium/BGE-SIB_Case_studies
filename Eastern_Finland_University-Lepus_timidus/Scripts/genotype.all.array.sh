#!/bin/bash

#SBATCH --job-name=hyr_arr_gen
#SBATCH --account=project_2002674
#SBATCH --time=72:00:00
#SBATCH --mem-per-cpu=20G
#SBATCH --cpus-per-task=2
#SBATCH --partition=small
#SBATCH --output=HyR_LTLEgt_arr_out_%j.txt
#SBATCH --error=HyR_LTLEgt_arr_err_%j.txt
#SBATCH --mail-type=END
#SBATCH --ntasks=1

###Array setup here
#SBATCH --array=1-44%5
#SBATCH --open-mode=truncate
#SBATCH --output=hybrid_combine_%a_output_%j.txt
#SBATCH --error=hybrid_combine_%a_errors_%j.txt

#modules
export PATH="/projappl/project_2002674/chromcomp_env//bin:$PATH"

#variables
refdir=/scratch/project_2002674/zsofia/Hybrid_rerun/GENOMES/
indir=/scratch/project_2002674/zsofia/Hybrid_rerun/04_Align/
vardir=/scratch/project_2002674/zsofia/Hybrid_rerun/05_Varcall/
gendir=/scratch/project_2002674/zsofia/Hybrid_rerun/06_Genotype/

chr=$(sed -n "${SLURM_ARRAY_TASK_ID}p" < chr_list.txt)
echo $chr

#LT
ref=$refdir/GCF_033115175.1_mLepTim1.pri_genomic.fna

gatk --java-options '-Xmx40G' CombineGVCFs -R $ref \
     --arguments_file all.vcf.list \
     -L $chr -O $gendir/combined.$chr.LepEur.g.vcf.gz

gatk --java-options '-Xmx40G' GenotypeGVCFs \
     -R $ref \
      -V $gendir/combined.$chr.LepEur.g.vcf.gz \
      -O $gendir/genotyped.$chr.LepEur.vcf.gz
