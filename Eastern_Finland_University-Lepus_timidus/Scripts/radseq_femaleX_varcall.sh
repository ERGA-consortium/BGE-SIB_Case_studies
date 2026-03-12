#!/bin/bash

#set up for small dataset
#SBATCH --job-name=radvarcall
#SBATCH --account=project_2002674
#SBATCH --time=72:00:00
#SBATCH --mem-per-cpu=10G
#SBATCH --cpus-per-task=8
#SBATCH --partition=small
#SBATCH --output=Radseq_varcall_output_%j.txt
#SBATCH --error=Radseq_varcall_errors_%j.txt
#SBATCH --mail-type=END
##SBATCH --gres=nvme:120

#modules
module load samtools
module load gatk

#variables
knownsites="lepeur_hoax.vcf"

workdir="/scratch/project_2002674/zsofia/Hybrid_hares/"

ref1="${workdir}/REFERENCE/GCF_033115175.1_mLepTim1.pri_genomic.fna"

out2="${workdir}/11_radseq_bqsr/"
vcfdir="${workdir}/17_haploid_calls/"
vault="/scratch/project_2002674/zsofia/Hybrid_hares/04_bwa_align/hybrids/"
logdir="${workdir}/script/log/"

#Make dirs
[ -d ${workdir}/17_haploid_calls ] || mkdir -p ${workdir}/17_haploid_calls

while read ID; do
    echo "Doing variant calling for ${ID} chr X..."
    gatk --java-options "-Xms70G -Xmx70G" HaplotypeCaller \
	 -ERC GVCF \
	 -I $out2/${ID}.bqsr.bam \
	 -O $vcfdir/${ID}.X.g.vcf \
	 -ploidy 2 \
	 -L NC_084850.1 \
	 -R $ref1 2> $logdir/${ID}.X.g.vcf.log
    echo "End work on $ID"
done < female.list


echo "Finished everything! \n Go check the output!"
