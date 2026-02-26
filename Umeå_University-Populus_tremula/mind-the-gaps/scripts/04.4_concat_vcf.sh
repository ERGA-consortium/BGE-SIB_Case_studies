#!/bin/bash -l

#SBATCH -A hpc2n2025-120	
#SBATCH -n 1
#SBATCH -c 4
#SBATCH -t 1-00:00:00
#SBATCH -J concat_vcf
#SBATCH --output=/proj/nobackup/hpc2nstor2025-059/mimmi/mind_the_gaps/reports/sbatch/concat_vcf/sbatch_R-%x_%j-%a.out
#SBATCH --error=/proj/nobackup/hpc2nstor2025-059/mimmi/mind_the_gaps/reports/sbatch/concat_vcf/sbatch_R-%x_%j-%a.err
#SBATCH -a 0-1

eval "$(mamba shell hook --shell bash)"
mamba activate chromcomp

cd /proj/nobackup/hpc2nstor2025-059/mimmi/mind_the_gaps

species=$1

run_list=(draft chrom)
asm=$(echo ${run_list[$SLURM_ARRAY_TASK_ID]})

echo ${species}
echo ${asm}

# list files
ls data/${species}/gvcf/${asm}/genomicDB/cohort_${asm}_*_raw.vcf.gz | \
    sort -V > data/${species}/gvcf/${asm}/genomicDB/list_vcf_to_concat.txt

# concat vcfs
bcftools concat \
    --threads 4 \
    --file-list data/${species}/gvcf/${asm}/genomicDB/list_vcf_to_concat.txt \
    --write-index \
    -Oz -o data/${species}/gvcf/${asm}/cohort.${asm}.raw.vcf.gz