#!/bin/bash
#SBATCH -c 16
#SBATCH --mem-per-cpu 2G
#SBATCH --time=8:00:00
#SBATCH --array=1-53%10
#SBATCH --job-name Rohan

dir=PATH

ref="PATH"

sample=$(sed -n "$SLURM_ARRAY_TASK_ID"p /home/qvw641/WhoopingCrane/list_bam_all.txt | awk '{print $1}')
name=$(sed -n "$SLURM_ARRAY_TASK_ID"p /home/qvw641/WhoopingCrane/Samples_all | awk '{print $1}')

# Historical
estimateDamage.pl --length 50 --threads 16 -o $dir/${name} $ref ${sample}
rohan -t 16  --size 500000 --rohmu 2e-5 --deam5p ${dir}/${name}.5p.prof --deam3p ${dir}/${name}.3p.prof --auto BigScaffolds_rohan.txt -o ${dir}/${name}_aDNA_2e5_500Kb_nochrZ $ref ${sample}

# modern without damage
rohan -t 16  --size 1000000 --rohmu 2e-5 --auto BigScaffolds_rohan.txt -o ${dir}/${name}_modern_2e5_1Mb_nochrZ $ref ${sample}