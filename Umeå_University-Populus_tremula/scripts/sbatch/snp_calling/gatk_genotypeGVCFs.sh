#!/bin/bash -l

#SBATCH -A hpc2n	
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 20
#SBATCH -t 1-00:00:00
#SBATCH -J genotypeGVCF
#SBATCH --output=reports/sbatch/genotypeGVCF/sbatch_R-%x_%j-%a.out
#SBATCH --error=reports/sbatch/genotypeGVCF/sbatch_R-%x_%j-%a.err
#SBATCH -a 19

ml GCCcore/12.3.0 
ml GATK/4.5.0.0-Java-17

chrm_ids=($(cat info_files/ids_chromosoms.txt))
chrm=$(echo ${chrm_ids[$SLURM_ARRAY_TASK_ID]})

ref="data/ref_genome/aspen/v2.2/fasta/Potra02_genome.fasta"

echo "gatk --java-options "-Xmx90g" GenotypeGVCFs -R ${ref} -V gendb://${chrm} -O ../chromosomes_vcfs/${chrm}.vcf.gz"
gatk --java-options "-Xmx90g" GenotypeGVCFs \
    -R ${ref} \
    -V gendb://${chrm} \
    -O ../chromosomes_vcfs/${chrm}.vcf.gz