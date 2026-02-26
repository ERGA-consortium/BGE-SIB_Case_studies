#!/bin/bash -l

#SBATCH -A hpc2n
#SBATCH -N 1	
#SBATCH -n 1
#SBATCH -c 1
#SBATCH -t 1-00:00:00
#SBATCH -J snp_calling_chrm
#SBATCH --output=reports/sbatch/snp_calling_chrm/sbatch_R-%x_%j-%a.out
#SBATCH --error=reports/sbatch/snp_calling_chrm/sbatch_R-%x_%j-%a.err
#SBATCH -a 0-18

ml GCCcore/12.3.0 
ml GATK/4.5.0.0-Java-17

chrm_ids=($(cat info_files/ids_chromosoms.txt))

sample=$1
chrm=$(echo ${chrm_ids[$SLURM_ARRAY_TASK_ID]})

data_dir="data/mapped/aspen/NordAsp/"
output_dir="data/variation/aspen/per_sample_runs/"${sample}/

ref="data/ref_genome/aspen/v2.2/fasta/Potra02_genome.fasta"

# if no sample dir, make one
if [ ! -d ${output_dir} ]; then
    mkdir -p ${output_dir}
fi

# step 1, haplotype caller
echo "gatk HaplotypeCaller -I ${data_dir}${sample}_sorted_RG_dedup.bam -O ${output_dir}${chrm}.g.vcf.gz -R ${ref} --emit-ref-confidence GVCF -L ${chrm} --native-pair-hmm-threads 1"
gatk HaplotypeCaller \
    -I ${data_dir}${sample}_sorted_RG_dedup.bam \
    -O ${output_dir}${chrm}.g.vcf.gz \
    -R ${ref} \
    --emit-ref-confidence GVCF \
    -L ${chrm} \
    --native-pair-hmm-threads 1