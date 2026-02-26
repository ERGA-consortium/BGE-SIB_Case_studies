#!/bin/bash -l

#SBATCH -A hpc2n2025-120	
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 20
#SBATCH -t 1-00:00:00
#SBATCH --mail-user mimmi.eriksson@slu.se
#SBATCH --mail-type=FAIL,END
#SBATCH -J genotypeGVCF
#SBATCH --output=/proj/nobackup/hpc2nstor2025-059/mimmi/aspen_snp_call/reports/sbatch/genotypeGVCF/sbatch_R-%x_%j-%a.out
#SBATCH --error=/proj/nobackup/hpc2nstor2025-059/mimmi/aspen_snp_call/reports/sbatch/genotypeGVCF/sbatch_R-%x_%j-%a.err
#SBATCH -a 19

ml GCCcore/12.3.0 
ml GATK/4.5.0.0-Java-17

cd /proj/nobackup/hpc2nstor2025-059/data/snp_calls/chromosomes_DBs

chrm_ids=($(cat /proj/nobackup/hpc2nstor2025-059/mimmi/aspen_snp_call/info_files/ids_chromosoms.txt))
chrm=$(echo ${chrm_ids[$SLURM_ARRAY_TASK_ID]})

ref="/proj/nobackup/hpc2nstor2025-059/data/ref_genome/aspen/v2.2/fasta/Potra02_genome.fasta"

echo "gatk --java-options "-Xmx90g" GenotypeGVCFs -R ${ref} -V gendb://${chrm} -O ../chromosomes_vcfs/${chrm}.vcf.gz"
gatk --java-options "-Xmx90g" GenotypeGVCFs \
    -R ${ref} \
    -V gendb://${chrm} \
    -O ../chromosomes_vcfs/${chrm}.vcf.gz