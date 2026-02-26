#!/bin/bash -l

#SBATCH -A hpc2n2025-120	
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 11
#SBATCH -t 1-00:00:00
#SBATCH --mail-user mimmi.eriksson@slu.se
#SBATCH --mail-type=FAIL,END
#SBATCH -J mapping
#SBATCH --output=/proj/nobackup/hpc2nstor2025-059/mimmi/aspen_snp_call/reports/sbatch/mapping/sbatch_R-%x_%j-%a.out
#SBATCH --error=/proj/nobackup/hpc2nstor2025-059/mimmi/aspen_snp_call/reports/sbatch/mapping/sbatch_R-%x_%j-%a.err
#SBATCH -a 0-39

ml GCC/13.2.0
ml BWA/0.7.18
ml SAMtools/1.19.2
ml picard/3.3.0-Java-17

cd /proj/nobackup/hpc2nstor2025-059

file_info=($(cat mimmi/aspen_snp_call/info_files/id_mapping_files.txt))

data_dir="data/fastq/aspen/NordAsp/trimmed/"
output_dir="data/mapped/aspen/NordAsp/"

ref="data/ref_genome/aspen/v2.2/fasta/Potra02_genome.fasta"

sample_id=$(echo ${file_info[$SLURM_ARRAY_TASK_ID]} | cut -d"_" -f1)
pu_info=$(echo ${file_info[$SLURM_ARRAY_TASK_ID]} | cut -d"_" -f2)

ID="NordAsp_H204SC25041259" # project
PU=${pu_info} # Format: flowcell.lane.library

SM="NordAsp_"${sample_id} # accession
LB="H204SC25041259_"${sample_id} # sample_id
PM="NovaSeq6000" # platform
PI=269 # Estimated insert size, output from fastp

outfile="NordAsp_"${sample_id}

r1=${data_dir}${sample_id}.R1.fq.gz
r2=${data_dir}${sample_id}.R2.fq.gz

# step 1, mapping
echo "bwa mem -t 10 -M ${ref} ${r1} ${r2} > ${output_dir}${outfile}.sam"
bwa mem -t 10 -M ${ref} ${r1} ${r2} > ${output_dir}${outfile}.sam

# step 2, flagstat
echo "samtools flagstat -O tsv ${output_dir}${outfile}.sam > mimmi/aspen_snp_call/reports/samtools_flagstats/${outfile}.tsv"
samtools flagstat -O tsv ${output_dir}${outfile}.sam > mimmi/aspen_snp_call/reports/samtools_flagstats/${outfile}.tsv

# step 3, samtools sort
echo "samtools view -b -u -F 4 -q 20 -@ 10 ${output_dir}${outfile}.sam | samtools sort -@ 10 - > ${output_dir}${outfile}_sorted.bam"
samtools view -b -u -F 4 -q 20 -@ 10 ${output_dir}${outfile}.sam | samtools sort -@ 10 - > ${output_dir}${outfile}_sorted.bam
samtools index ${output_dir}${outfile}_sorted.bam
rm ${output_dir}${outfile}.sam

# step 4 change @RG info
echo "java -jar $EBROOTPICARD/picard.jar AddOrReplaceReadGroups I=${output_dir}${outfile}_sorted.bam O=${output_dir}${outfile}_sorted_RG.bam RGID=${ID} RGLB=${LB} RGPL=ILLUMINA RGPU=${PU} RGSM=${SM} RGPM=${PM} RGPI=${PI} CREATE_INDEX=true"
java -jar $EBROOTPICARD/picard.jar AddOrReplaceReadGroups \
    I=${output_dir}${outfile}_sorted.bam \
    O=${output_dir}${outfile}_sorted_RG.bam \
    RGID=${ID} \
    RGLB=${LB} \
    RGPL=ILLUMINA \
    RGPU=${PU} \
    RGSM=${SM} \
    RGPM=${PM} \
    RGPI=${PI} \
    CREATE_INDEX=true

rm ${output_dir}${outfile}_sorted.bam
rm ${output_dir}${outfile}_sorted.bai

# step 5, picard alignment summary metrics
echo "java -jar $EBROOTPICARD/picard.jar CollectAlignmentSummaryMetrics R=$ref I=${output_dir}${outfile}_sorted_RG.bam O=mimmi/aspen_snp_call/reports/picard_alignment_summary/${outfile}_alignment_metrics.txt"
java -jar $EBROOTPICARD/picard.jar CollectAlignmentSummaryMetrics \
    R=$ref \
    I=${output_dir}${outfile}_sorted_RG.bam \
    O=mimmi/aspen_snp_call/reports/picard_alignment_summary/${outfile}_alignment_metrics.txt

# step 6, picard mark duplicates
java -Xmx48g -jar $EBROOTPICARD/picard.jar MarkDuplicates \
    INPUT=${output_dir}${outfile}_sorted_RG.bam \
    OUTPUT=${output_dir}${outfile}_sorted_RG_dedup.bam \
    METRICS_FILE=mimmi/aspen_snp_call/reports/picard_mark_duplicates_metrics/${outfile}.picard_dup_metrics.txt \
    CREATE_INDEX=true

rm ${output_dir}${outfile}_sorted_RG.bam
rm ${output_dir}${outfile}_sorted_RG.bai