#!/bin/bash
#SBATCH -c 24
#SBATCH --mem-per-cpu 18G
#SBATCH --time=30:00:00
#SBATCH --job-name Het_1st

module load angsd/0.940

ref="/projects/echo/data/OWL/data/refGenome_NEW/bAthNoc2.reordered.hap1.fasta"

dir=${1}
bams=${2}
scaffolds=${3}
minind=${4}
trans=${5}
minMaf=${6}


angsd -b $bams -ref $ref -rf $scaffolds -anc $ref  -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -rmTrans $trans -trim 0 -C 50 -baq 0 -minMaf $minMaf -minMapQ 20 -minQ 20 -skipTriallelic 1 -doCounts 1 -doMajorMinor 1 -GL 2 -doGlf 2 -doMaf 2 -minInd $minind -out ${dir}Species_pop_filter_minind${minind}_baq${baq}

zcat ${dir}Species_pop_filter_minind${minind}_baq${baq}.mafs.gz  | awk '{print $1 , $2}' | sed 1d > ${dir}Species_pop_filter_minind${minind}_baq${baq}.list

angsd sites index ${dir}Species_pop_filter_minind${minind}_baq${baq}.list