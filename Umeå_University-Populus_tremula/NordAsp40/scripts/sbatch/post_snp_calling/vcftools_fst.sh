#!/bin/bash -l

#SBATCH -A hpc2n2025-120	
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 1
#SBATCH -t 1-00:00:00
#SBATCH -J vcftools_fst
#SBATCH --output=/home/m/mimmie/projects/aspen_hpc2nstor2025-059/mimmi/NordAsp_40/reports/sbatch/vcftools_fst/sbatch_R-%x_%j.out
#SBATCH --error=/home/m/mimmie/projects/aspen_hpc2nstor2025-059/mimmi/NordAsp_40/reports/sbatch/vcftools_fst/sbatch_R-%x_%j.err

ml GCC/13.2.0
ml VCFtools/0.1.16

cd /proj/nobackup/hpc2nstor2025-059/mimmi/NordAsp_40

input_vcf="output_data/vcf_filtered/NordAsp_H204SC25041259_biallelic-missing0.8-maf0.05_recoded_LDpruned.vcf.gz"

# Pairwise FST (two example populations)
vcftools --gzvcf ${input_vcf} \
     	--weir-fst-pop output_data/pop_files/east.txt \
     	--weir-fst-pop output_data/pop_files/mountain.txt \
      --fst-window-size 10000 \
     	--out output_data/vcftools/NordAsp40_Fst_east-mountain

# Pairwise FST (two example populations)
vcftools --gzvcf ${input_vcf} \
     	--weir-fst-pop output_data/pop_files/east.txt \
      --weir-fst-pop output_data/pop_files/west.txt \
      --fst-window-size 10000 \
     	--out output_data/vcftools/NordAsp40_Fst_east-west

# Pairwise FST (two example populations)
vcftools --gzvcf ${input_vcf} \
     	--weir-fst-pop output_data/pop_files/mountain.txt \
      --weir-fst-pop output_data/pop_files/west.txt \
      --fst-window-size 10000 \
     	--out output_data/vcftools/NordAsp40_Fst_mountain-west