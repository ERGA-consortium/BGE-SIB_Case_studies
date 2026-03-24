#!/bin/bash
#SBATCH -c 10
#SBATCH --mem-per-cpu 8G
#SBATCH --time=12:00:00
#SBATCH --partition=cpuqueue
#SBATCH --job-name pca-ldt126

module load angsd/0.940
module load pcangsd/1.36.3

ref="/projects/echo/data/OWL/data/refGenome_NEW/bAthNoc2.reordered.hap1.fasta"


dir=${1}
scaffolds=${2}
bams=${3}
minInd=${4}
minMaf=${5}
trans=${6}

angsd -bam $bams -ref $ref -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -doCounts 1 -rmTrans $trans -C 50 -baq 0 -minInd 10 -skipTriallelic 1 -GL 2 -minQ 20 -minMapQ 30 -nThreads 10 -doGlf 2 -doMajorMinor 1 -doMaf 1 -doHWE 1 -SNP_pval 1e-6 -minMaf $minMaf -dosnpstat 1 -HWE_pval 1e-2 -rf $scaffolds -out ${dir}pca_out_minind${minInd}_trans${trans}_minMaf${minMaf}

pcangsd -b ${dir}*.beagle.gz -o ${dir}out_pca_Ind${minInd}_trans${trans}.pcangsd
