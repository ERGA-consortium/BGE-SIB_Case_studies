#!/bin/bash

#SBATCH --job-name=Ycombine
#SBATCH --account=project_2002674
#SBATCH --time=72:00:00
#SBATCH --partition=small
#SBATCH --mail-type=END

#SBATCH --mem-per-cpu=20G
#SBATCH --cpus-per-task=2
#SBATCH --output=Y_combine_output_%j.txt
#SBATCH --error=Y_combine_errors_%j.txt

#modules
module load samtools
module load gatk

#variables
knownsites="lepeur_hoax.vcf"
workdir="/scratch/project_2002674/zsofia/Hybrid_hares/"
ref1="${workdir}/REFERENCE/GCF_033115175.1_mLepTim1.pri_genomic.fna"

logdir="${workdir}/script/log/"
outdir="${workdir}/18_combined_xymt/"
chr="NC_084851.1"


[ -d $outdir ] || mkdir -p $outdir


#Combine and filter GVCFs
echo "Running CombineGVCFs for ${chr}..."
srun gatk --java-options "-Xmx18G -Xms18G" CombineGVCFs \
     -O $outdir/Y.combined.g.vcf.gz \
     -L $chr \
     -R $ref1 \
     --arguments_file ycombine2.args.list \
     2> $logdir/2024-08-Y-combinegvcfs.log

exit
