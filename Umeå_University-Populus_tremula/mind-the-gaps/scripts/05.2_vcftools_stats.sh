#!/bin/bash -l

#SBATCH -A hpc2n2025-120	
#SBATCH -n 1
#SBATCH -c 1
#SBATCH -t 1-00:00:00
#SBATCH -J vcftools_stats
#SBATCH --output=/proj/nobackup/hpc2nstor2025-059/mimmi/mind_the_gaps/reports/sbatch/vcftools_stats/sbatch_R-%x_%j-%a.out
#SBATCH --error=/proj/nobackup/hpc2nstor2025-059/mimmi/mind_the_gaps/reports/sbatch/vcftools_stats/sbatch_R-%x_%j-%a.err
#SBATCH -a 0-1

eval "$(mamba shell hook --shell bash)"
mamba activate chromcomp

cd /proj/nobackup/hpc2nstor2025-059/mimmi/mind_the_gaps

run_list=(draft chrom)
asm=$(echo ${run_list[$SLURM_ARRAY_TASK_ID]})

species=$1

echo ${species}
echo ${asm}

# Nucleotide diversity (π)
vcftools --gzvcf data/${species}/gvcf/${asm}/cohort.${asm}.final.recoded.vcf.gz \
     	--site-pi \
     	--out data/${species}/stats/${asm}/pi_${asm}

# Observed heterozygosity
vcftools --gzvcf data/${species}/gvcf/${asm}/cohort.${asm}.final.recoded.vcf.gz \
     	--het \
     	--out data/${species}/stats/${asm}/het_${asm}

# Pairwise FST (two example populations)
vcftools --gzvcf data/${species}/gvcf/${asm}/cohort.${asm}.final.recoded.vcf.gz \
     	--weir-fst-pop data/${species}/pop1.txt \
     	--weir-fst-pop data/${species}/pop2.txt \
     	--out data/${species}/stats/${asm}/fst_${asm}