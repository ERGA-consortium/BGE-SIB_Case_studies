#!/bin/bash -l

#SBATCH -A hpc2n2025-120	
#SBATCH -n 1
#SBATCH -c 4
#SBATCH -t 7-00:00:00
#SBATCH -J haplo_merge
#SBATCH --output=/proj/nobackup/hpc2nstor2025-059/mimmi/mind_the_gaps/reports/sbatch/haplo_merge/sbatch_R-%x_%j-%a.out
#SBATCH --error=/proj/nobackup/hpc2nstor2025-059/mimmi/mind_the_gaps/reports/sbatch/haplo_merge/sbatch_R-%x_%j-%a.err
#SBATCH -a 0-19 

eval "$(mamba shell hook --shell bash)"
mamba activate chromcomp

cd /proj/nobackup/hpc2nstor2025-059/mimmi/mind_the_gaps

species=$1
asm=$2

echo ${species}
echo ${asm}

bed_list=($(ls data/${species}/bedfiles/${asm}/*.bed))
bed_file=$(echo ${bed_list[$SLURM_ARRAY_TASK_ID]})

run_id=$(basename ${bed_file/.bed//} | cut -d"_" -f2)

## Make cohort file
# ex:
# sample1	/path/to/sample1.g.vcf.gz
# sample2	/path/to/sample2.g.vcf.gz
# sample3	/path/to/sample3.g.vcf.gz

# list all files and add the full path
ls data/${species}/gvcf/${asm}/sample_calls/*_${asm}_${run_id}.g.vcf.gz | \
    sort -V | \
    sed 's/^/\/proj\/nobackup\/hpc2nstor2025-059\/mimmi\/mind_the_gaps\//g' > data/${species}/gvcf/${asm}/genomicDB/${asm}_${run_id}_cohort-files.txt

# add sample id in the fist column and save into a new temp file
cat data/populus_tremula/sample_info.txt | \
    cut -f1 | \
    sort -V | \
    paste - data/${species}/gvcf/${asm}/genomicDB/${asm}_${run_id}_cohort-files.txt > data/${species}/gvcf/${asm}/genomicDB/${asm}_${run_id}_cohort-files_temp.txt

# move the temp file to the old file name
mv data/${species}/gvcf/${asm}/genomicDB/${asm}_${run_id}_cohort-files_temp.txt data/${species}/gvcf/${asm}/genomicDB/${asm}_${run_id}_cohort-files.txt

## genomicsDBimport, combine sample calls
gatk --java-options "-Xmx20g -DGATK_STACKTRACE_ON_USER_EXCEPTION=true" GenomicsDBImport \
    --genomicsdb-workspace-path data/${species}/gvcf/${asm}/genomicDB/range_${run_id} \
    --sample-name-map data/${species}/gvcf/${asm}/genomicDB/${asm}_${run_id}_cohort-files.txt \
	--merge-contigs-into-num-partitions 25 \
    --merge-input-intervals \
	--bypass-feature-reader \
	-L ${bed_file}