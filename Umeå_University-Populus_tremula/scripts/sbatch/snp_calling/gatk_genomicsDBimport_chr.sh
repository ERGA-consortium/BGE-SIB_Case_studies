#!/bin/bash -l

#SBATCH -A hpc2n	
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 25
#SBATCH -t 1-00:00:00
#SBATCH -J genomicDB_chr
#SBATCH --output=reports/sbatch/genomicDB_chr/sbatch_R-%x_%j-%a.out
#SBATCH --error=reports/sbatch/genomicDB_chr/sbatch_R-%x_%j-%a.err
#SBATCH -a 0-18

ml GCCcore/12.3.0 
ml GATK/4.5.0.0-Java-17

chrm_ids=($(cat info_files/ids_chromosoms.txt))

chrm=$(echo ${chrm_ids[$SLURM_ARRAY_TASK_ID]})

data_dir="data/variation/aspen/"
output_dir=${data_dir}"chromosomes_DBs/"${chrm}

ref="data/ref_genome/aspen/v2.2/fasta/Potra02_genome.fasta"

## Make cohort file
# ex:
# NordAsp_A26	data/variation/aspen/per_sample_runs/NordAsp_A26/chr1.g.vcf.gz
# NordAsp_B7	data/variation/aspen/per_sample_runs/NordAsp_B7/chr1.g.vcf.gz
# NordAsp_C42	data/variation/aspen/per_sample_runs/NordAsp_C42/chr1.g.vcf.gz
# NordAsp_C46	data/variation/aspen/per_sample_runs/NordAsp_C46/chr1.g.vcf.gz
# NordAsp_D4	data/variation/aspen/per_sample_runs/NordAsp_D4/chr1.g.vcf.gz

# list all files and add the full path
ls ${data_dir}per_sample_runs/*/${chrm}.g.vcf.gz | \
    sed 's/^/\/proj\/nobackup\/hpc2nstor2025-059\//' > info_files/cohort_files/${chrm}.txt

# add sample id in the fist column and save into a new temp file
cat info_files/cohort_files/${chrm}.txt | \
    awk '{sub(/.*_runs\//, ""); sub(/\/.*/, ""); print}' | \
    paste - info_files/cohort_files/${chrm}.txt > info_files/cohort_files/${chrm}_temp.txt

# move the temp file to the old file name
mv info_files/cohort_files/${chrm}_temp.txt info_files/cohort_files/${chrm}.txt

## genomicsDBimport
echo "gatk --java-options "-Xmx120g -DGATK_STACKTRACE_ON_USER_EXCEPTION=true" GenomicsDBImport --genomicsdb-workspace-path ${output_dir} --sample-name-map info_files/cohort_files/${chrm}.txt --batch-size 20 --merge-contigs-into-num-partitions 25 --bypass-feature-reader -L ${chrm}"
gatk --java-options "-Xmx120g -DGATK_STACKTRACE_ON_USER_EXCEPTION=true" GenomicsDBImport \
    --genomicsdb-workspace-path ${output_dir} \
    --sample-name-map info_files/cohort_files/${chrm}.txt \
	--batch-size 20 \
	--merge-contigs-into-num-partitions 25 \
	--bypass-feature-reader \
	-L ${chrm}