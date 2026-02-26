#!/bin/bash -l

#SBATCH -A hpc2n2025-120	
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 25
#SBATCH -t 1-00:00:00
#SBATCH --mail-user mimmi.eriksson@slu.se
#SBATCH --mail-type=FAIL,END
#SBATCH -J genomicDB_scaffolds
#SBATCH --output=/proj/nobackup/hpc2nstor2025-059/mimmi/aspen_snp_call/reports/sbatch/genomicDB_scaffolds/sbatch_R-%x_%j.out
#SBATCH --error=/proj/nobackup/hpc2nstor2025-059/mimmi/aspen_snp_call/reports/sbatch/genomicDB_scaffolds/sbatch_R-%x_%j.err

ml GCCcore/12.3.0 
ml GATK/4.5.0.0-Java-17

cd /proj/nobackup/hpc2nstor2025-059

data_dir="data/variation/aspen/"
output_dir=${data_dir}"chromosomes_DBs/scaffolds"

ref="/proj/nobackup/hpc2nstor2025-059/data/ref_genome/aspen/v2.2/fasta/Potra02_genome.fasta"

## Make cohort file
# ex:
# NordAsp_A26	/proj/nobackup/hpc2nstor2025-059/data/variation/aspen/per_sample_runs/NordAsp_A26/chr1.g.vcf.gz
# NordAsp_B7	/proj/nobackup/hpc2nstor2025-059/data/variation/aspen/per_sample_runs/NordAsp_B7/chr1.g.vcf.gz
# NordAsp_C42	/proj/nobackup/hpc2nstor2025-059/data/variation/aspen/per_sample_runs/NordAsp_C42/chr1.g.vcf.gz
# NordAsp_C46	/proj/nobackup/hpc2nstor2025-059/data/variation/aspen/per_sample_runs/NordAsp_C46/chr1.g.vcf.gz
# NordAsp_D4	/proj/nobackup/hpc2nstor2025-059/data/variation/aspen/per_sample_runs/NordAsp_D4/chr1.g.vcf.gz

# list all files and add the full path
ls ${data_dir}per_sample_runs/*/scaffolds.g.vcf.gz | \
    sed 's/^/\/proj\/nobackup\/hpc2nstor2025-059\//' > mimmi/aspen_snp_call/info_files/cohort_files/scaffolds.txt

# add sample id in the fist column and save into a new temp file
cat mimmi/aspen_snp_call/info_files/cohort_files/scaffolds.txt | \
    awk '{sub(/.*_runs\//, ""); sub(/\/.*/, ""); print}' | \
    paste - mimmi/aspen_snp_call/info_files/cohort_files/scaffolds.txt > mimmi/aspen_snp_call/info_files/cohort_files/scaffolds_temp.txt

# move the temp file to the old file name
mv mimmi/aspen_snp_call/info_files/cohort_files/scaffolds_temp.txt mimmi/aspen_snp_call/info_files/cohort_files/scaffolds.txt

## genomicsDBimport
echo "gatk --java-options "-Xmx120g -DGATK_STACKTRACE_ON_USER_EXCEPTION=true" GenomicsDBImport --genomicsdb-workspace-path ${output_dir} --sample-name-map mimmi/aspen_snp_call/info_files/cohort_files/scaffolds.txt --batch-size 20 --merge-contigs-into-num-partitions 25 --bypass-feature-reader -L /proj/nobackup/hpc2nstor2025-059/mimmi/aspen_snp_call/info_files/scaffolds.bed"
gatk --java-options "-Xmx120g -DGATK_STACKTRACE_ON_USER_EXCEPTION=true" GenomicsDBImport \
    --genomicsdb-workspace-path ${output_dir} \
    --sample-name-map mimmi/aspen_snp_call/info_files/cohort_files/scaffolds.txt \
	--batch-size 20 \
	--merge-contigs-into-num-partitions 25 \
	--bypass-feature-reader \
	-L /proj/nobackup/hpc2nstor2025-059/mimmi/aspen_snp_call/info_files/scaffolds.bed