#!/bin/bash -l

#SBATCH -A hpc2n2025-120	
#SBATCH -n 1
#SBATCH -c 4
#SBATCH -t 4-00:00:00
#SBATCH -J genotype_vcf
#SBATCH --output=/proj/nobackup/hpc2nstor2025-059/mimmi/mind_the_gaps/reports/sbatch/genotype_vcf/sbatch_R-%x_%j-%a.out
#SBATCH --error=/proj/nobackup/hpc2nstor2025-059/mimmi/mind_the_gaps/reports/sbatch/genotype_vcf/sbatch_R-%x_%j-%a.err
#SBATCH -a 0-19

eval "$(mamba shell hook --shell bash)"
mamba activate chromcomp

cd /proj/nobackup/hpc2nstor2025-059/mimmi/mind_the_gaps

species=$1
asm=$2

echo ${species}
echo ${asm}

bed_list=($(ls data/${species}/bedfiles/${asm}/*.bed))
bed_file=$(echo ${bed_list[$SLURM_ARRAY_TASK_ID]})

run_id=$(basename ${bed_file/.bed//} | cut -d"_" -f2)

ref="/proj/nobackup/hpc2nstor2025-059/mimmi/mind_the_gaps/data/${species}/assemblies/${asm}/${asm}_round-2_masked.fa"

cd data/${species}/gvcf/${asm}/genomicDB

# Cohort‑level genotyping
gatk --java-options "-Xmx20g" GenotypeGVCFs \
    -R ${ref} \
    -V gendb://range_${run_id} \
    -O cohort_${asm}_${run_id}_raw.vcf.gz