#!/bin/bash -l

#SBATCH -A hpc2n2025-120
#SBATCH -N 1	
#SBATCH -n 1
#SBATCH -c 2
#SBATCH -t 1-00:00:00
#SBATCH --mail-user mimmi.eriksson@slu.se
#SBATCH --mail-type=FAIL,END
#SBATCH -J snp_calling_scaffolds
#SBATCH --output=/proj/nobackup/hpc2nstor2025-059/mimmi/aspen_snp_call/reports/sbatch/snp_calling_scaffolds/sbatch_R-%x_%j-%a.out
#SBATCH --error=/proj/nobackup/hpc2nstor2025-059/mimmi/aspen_snp_call/reports/sbatch/snp_calling_scaffolds/sbatch_R-%x_%j-%a.err
#SBATCH -a 0-39

ml GCCcore/12.3.0 
ml GATK/4.5.0.0-Java-17

cd /proj/nobackup/hpc2nstor2025-059

sample_ids=($(cat mimmi/aspen_snp_call/info_files/sample_names.txt))
sample=$(echo ${sample_ids[$SLURM_ARRAY_TASK_ID]})

data_dir="data/mapped/aspen/NordAsp/"
output_dir="data/variation/aspen/per_sample_runs/"${sample}/

ref="data/ref_genome/aspen/v2.2/fasta/Potra02_genome.fasta"

# if no sample dir, make one
if [ ! -d ${output_dir} ]; then
    mkdir -p ${output_dir}
fi

# step 1, haplotype caller
echo "gatk HaplotypeCaller -I ${data_dir}${sample}_sorted_RG_dedup.bam -O ${output_dir}scaffolds.g.vcf.gz -R ${ref} --emit-ref-confidence GVCF -L mimmi/aspen_snp_call/info_files/scaffolds.bed --native-pair-hmm-threads 1"
gatk HaplotypeCaller \
    -I ${data_dir}${sample}_sorted_RG_dedup.bam \
    -O ${output_dir}scaffolds.g.vcf.gz \
    -R ${ref} \
    --emit-ref-confidence GVCF \
    -L mimmi/aspen_snp_call/info_files/scaffolds.bed \
    --native-pair-hmm-threads 1