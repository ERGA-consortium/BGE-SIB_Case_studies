#!/bin/bash
#SBATCH -c 10
#SBATCH --mem-per-cpu 2G
#SBATCH --time=4:00:00
#SBATCH --array=1-25%10
#SBATCH --job-name Het_2nd

module load angsd/0.940

dir=${1}
ref=${2}
bams=${3}
baq=${4}
mindepth=${5}
trans=${6}

sample=$(sed -n "$SLURM_ARRAY_TASK_ID"p $bams | awk '{print $1}')
# extract basename and remove trailing
name=$(sed -n "$SLURM_ARRAY_TASK_ID"p $bams | awk '{print $1}' | xargs -n1 basename -s .rmdup.bam)


angsd  -i ${sample} -ref $ref -sites ${dir}Species_pop_filter*.list -anc $ref -out ${dir}${name}_filtList_pop -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -noTrans $trans -C 50 -baq $baq -minMapQ 30 -minQ 20  -setMinDepth $mindepth -setMaxDepth 50 -doCounts 1 -nThreads 10 -GL 2 -doSaf 1 -rf /projects/echo/people/ldt126/dirs/scafflist_100Kb.txt
#try -baq with 1 and 2

realSFS -fold 1 -bootstrap 5 ${dir}${name}_filtList_pop.saf.idx > ${dir}out_het_${name}_minDpthIndv${mindepth}_trans${trans}_baq${baq}_angsd.est.ml
