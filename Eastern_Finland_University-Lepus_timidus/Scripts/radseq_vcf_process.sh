#!/bin/bash

#set up for small dataset
#SBATCH --job-name=radvcf
#SBATCH --account=project_2002674
#SBATCH --time=72:00:00
#SBATCH --mem-per-cpu=10G
#SBATCH --cpus-per-task=4
#SBATCH --partition=small
#SBATCH --output=hybrid_filter_output_%j.txt
#SBATCH --error=hybrid_filter_errors_%j.txt
#SBATCH --mail-type=END

#modules
module load samtools/1.21
module load gatk/4.5.0.0

#variables
knownsites="lepeur_hoax.vcf"

workdir="/scratch/project_2002674/zsofia/Hybrid_hares/"

ref1="${workdir}/REFERENCE/GCF_033115175.1_mLepTim1.pri_genomic.fna"

logdir="${workdir}/script/log/"
outdir="${workdir}/06_Genotype/"


#Filtering combined and genotyped vcf
echo "Separate variants for filtration..."
gatk --java-options "-Xmx40G" SelectVariants \
     -V $outdir/radseq.genotyped.vcf.gz \
     -select-type SNP \
     -O $outdir/radseq.snps.vcf.gz \
     2> $logdir/selectsnps.log
gatk --java-options "-Xmx110G" SelectVariants \
     -V $outdir/radseq.genotyped.vcf.gz \
     -select-type INDEL \
     -select-type MIXED \
     -O $outdir/radseq.indel.vcf.gz \
     2> $logdir/selectindels.log

echo "RUnning VariantFiltration..."
gatk --java-options "-Xmx40G" VariantFiltration \
     -V $outdir/radseq.snps.vcf.gz \
    -filter "QD < 2.0" --filter-name "QD2" \
    -filter "QUAL < 30.0" --filter-name "QUAL30" \
    -filter "SOR > 3.0" --filter-name "SOR3" \
    -filter "FS > 60.0" --filter-name "FS60" \
    -filter "MQ < 40.0" --filter-name "MQ40" \
    -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
    -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
    -O $outdir/radseq.snps.filtered.vcf.gz \
    2> $logdir/filtersnps.log

gatk --java-options "-Xmx40G" VariantFiltration \
     -V $outdir/radseq.indel.vcf.gz \
     -filter "QD < 2.0" --filter-name "QD2" \
     -filter "QUAL < 30.0" --filter-name "QUAL30" \
     -filter "FS > 200.0" --filter-name "FS200" \
     -filter "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum-20" \
     -O $outdir/radseq.indels.filtered.vcf.gz \
     2> $logdir/filterindels.log

echo "Merging outputs..."
gatk --java-options "-Xmx40G" MergeVcfs \
     -I $outdir/radseq.snps.filtered.vcf.gz \
     -I $outdir/radseq.indels.filtered.vcf.gz \
     -O $outdir/radseq.all.filtered.vcf.gz
