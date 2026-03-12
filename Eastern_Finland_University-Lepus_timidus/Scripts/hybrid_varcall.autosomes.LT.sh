#!/bin/bash

#SBATCH --job-name=hyvarcall
#SBATCH --account=project_2002674
#SBATCH --time=72:00:00
#SBATCH --mem-per-cpu=4G
#SBATCH --cpus-per-task=3
#SBATCH --partition=test
#SBATCH --output=hybrid_varcall_out_%j.txt
#SBATCH --error=hybrid_varcall_err_%j.txt
#SBATCH --mail-type=END

###Array setup here
##SBATCH --array=1-20%5
##SBATCH --open-mode=truncate

#modules
module load gatk/4.5.0.0

#variables
refdir=/scratch/project_2002674/zsofia/Hybrid_hares/REFERENCE/
indir=/scratch/project_2002674/zsofia/Hybrid_hares/04_Align/
vardir=/scratch/project_2002674/zsofia/Hybrid_hares/05_Varcall/
gendir=/scratch/project_2002674/zsofia/Hybrid_hares/06_Genotype/

#LE
ref=$refdir/GCF_033115175.1_mLepTim1.pri_genomic.fna

chr=$(sed -n "${SLURM_ARRAY_TASK_ID}p" < chr_list.txt)

while read sample_id; do
    if [ ! -f $vardir/${sample_id}.LepEur.g.vcf.gz.tbi ]; then
	bamf=$indir/${sample_id}.LepEur.bam
	
	gatk --java-options '-Xmx8G' HaplotypeCaller \
     	     -R $ref \
     	     -I $bamf \
	     -L $chr \
     	     -ERC GVCF \
     	     -O $vardir/${sample_id}.LepEur.g.vcf.gz;
    else
	echo "$vardir/${sample_id}.LepEur.g.vcf.gz.tbi exists"
    fi
done < LT.confirmed.samples.list #Repeat with LT
