#!/bin/bash

#set up for small dataset
#SBATCH --job-name=hy_bwa
#SBATCH --account=project_2002674
#SBATCH --time=24:00:00
#SBATCH --mem-per-cpu=10G
#SBATCH --cpus-per-task=8
#SBATCH --partition=small
#SBATCH --output=hybrid_bwa_output_%j.txt
#SBATCH --error=hybrid_bwa_errors_%j.txt
#SBATCH --mail-type=END

#modules
module load samtools/1.21
module load bwa-mem2/2.2.1

#variables
workdir="/scratch/project_2002674/zsofia/Hybrid_hares/"

ref1="${workdir}/REFERENCE/GCF_033115175.1_mLepTim1.pri_genomic.fna"

out1="${workdir}/04_bwa_align/"
vault="${workdir}/03_cleanfq/"

#Make dirs
[ -d ${workdir}/04_bwa_align ] || mkdir ${workdir}/04_bwa_align

#Index
#echo "Creating indexes for ${ref1}"
#samtools faidx $ref1
#bwa-mem2 index $ref1

#Align
echo "Starting alingment..."

while read sample; do
    echo $sample
    infq1="${vault}/${sample}.1.fastq.gz"
    infq2="${vault}/${sample}.2.fastq.gz"
    bwa-mem2 mem -t 16 -M \
	     -R "@RG\tID:${sample}paired\tSM:${sample}\tPL:illumina" \
	     $ref1 $infq1 $infq2 | \
    samtools sort -m 8G -@ 4 | \
    samtools view -b -F 4 -m 8G -@ 4 > ${out1}/${sample}.bam ; \
    samtools index ${out1}/${sample}.bam;
done < samples.list

echo "Finished alignment!"
