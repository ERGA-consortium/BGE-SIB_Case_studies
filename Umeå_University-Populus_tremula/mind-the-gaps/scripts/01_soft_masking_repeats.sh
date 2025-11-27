#!/bin/bash -l

#SBATCH -A hpc2n2025-120
#SBATCH -N 1	
#SBATCH -n 1
#SBATCH -c 100
#SBATCH -t 7-00:00:00
#SBATCH -J masking_repeats
#SBATCH --output=/proj/nobackup/hpc2nstor2025-059/mimmi/mind_the_gaps/reports/sbatch/masking_repeats/sbatch_R-%x_%j-%a.out
#SBATCH --error=/proj/nobackup/hpc2nstor2025-059/mimmi/mind_the_gaps/reports/sbatch/masking_repeats/sbatch_R-%x_%j-%a.err
##SBATCH -a 0-1

eval "$(mamba shell hook --shell bash)"
mamba activate chromcomp

cd /proj/nobackup/hpc2nstor2025-059/mimmi/mind_the_gaps

species=$1
currated_lib=$2
# run_list=(draft chrom)
# asm=$(echo ${run_list[$SLURM_ARRAY_TASK_ID]})
asm="draft"

echo ${asm}
echo ${species}

### round one with the common plants database
# -pa parallel threads 
# -species "viridiplantae" repeat library hint (adjust to your clade)
# -xsmall  soft‑mask: convert repeats to lowercase
# -gff  emit annotation for downstream stats
# -dir output directory (writes *.masked.fa)
RepeatMasker \
	-pa 28 \
	-species "viridiplantae" \
	-xsmall \
	-gff \
	-dir data/${species}/assemblies/${asm}/ \
	data/${species}/assemblies/${asm}/${asm}.fa

### rename files produced in round one
mv data/${species}/assemblies/${asm}/${asm}.fa.out.gff data/${species}/assemblies/${asm}/${asm}_round-1.out.gff
mv data/${species}/assemblies/${asm}/${asm}.fa.cat.gz data/${species}/assemblies/${asm}/${asm}_round-1.cat.gz
mv data/${species}/assemblies/${asm}/${asm}.fa.masked data/${species}/assemblies/${asm}/${asm}_round-1_masked.fa
mv data/${species}/assemblies/${asm}/${asm}.fa.out data/${species}/assemblies/${asm}/${asm}_round-1.out
mv data/${species}/assemblies/${asm}/${asm}.fa.tbl data/${species}/assemblies/${asm}/${asm}_round-1.tbl

### round two with speces currated repeat library
RepeatMasker \
	-pa 28 \
	-engine ncbi \
	-lib ${currated_lib} \
	-xsmall \
	-gff \
	-dir data/${species}/assemblies/${asm}/ \
	data/${species}/assemblies/${asm}/${asm}_round-1_masked.fa

### rename files produced in round one
mv data/${species}/assemblies/${asm}/${asm}_round-1_masked.fa.out.gff data/${species}/assemblies/${asm}/${asm}_round-2.out.gff
mv data/${species}/assemblies/${asm}/${asm}_round-1_masked.fa.cat.gz data/${species}/assemblies/${asm}/${asm}_round-2.cat.gz
mv data/${species}/assemblies/${asm}/${asm}_round-1_masked.fa.masked data/${species}/assemblies/${asm}/${asm}_round-2_masked.fa
mv data/${species}/assemblies/${asm}/${asm}_round-1_masked.fa.out data/${species}/assemblies/${asm}/${asm}_round-2.out
mv data/${species}/assemblies/${asm}/${asm}_round-1_masked.fa.tbl data/${species}/assemblies/${asm}/${asm}_round-2.tbl
