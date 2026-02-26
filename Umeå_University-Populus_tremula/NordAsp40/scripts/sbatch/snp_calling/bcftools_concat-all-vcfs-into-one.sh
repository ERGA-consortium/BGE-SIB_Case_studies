#!/bin/bash -l

#SBATCH -A hpc2n2025-120	
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 10
#SBATCH -t 1-00:00:00
#SBATCH --mail-user mimmi.eriksson@slu.se
#SBATCH --mail-type=FAIL,END
#SBATCH -J concatVCFs
#SBATCH --output=/proj/nobackup/hpc2nstor2025-059/mimmi/aspen_snp_call/reports/sbatch/concatVCFs/sbatch_R-%x_%j.out
#SBATCH --error=/proj/nobackup/hpc2nstor2025-059/mimmi/aspen_snp_call/reports/sbatch/concatVCFs/sbatch_R-%x_%j.err

ml GCC/13.2.0 
ml BCFtools/1.19

cd /proj/nobackup/hpc2nstor2025-059/

bcftools concat \
    --threads 10 \
    --file-list /proj/nobackup/hpc2nstor2025-059/mimmi/aspen_snp_call/info_files/all_vcfs_to_concat.txt \
    --write-index \
    -Oz -o /proj/nobackup/hpc2nstor2025-059/data/variation/aspen/NordAsp_H204SC25041259.vcf.gz