#!/bin/bash

#Identify reads with primers at both ends and remove them
for i in *bam; do 
	lima "$i" primers.fa "${i/_reads.bam/.noPrimers_reads.bam}" --overwrite-biosample-names --isoseq --peek-guess;
done
#dentify and Remove polyA tails and artificial concatemers
for i in *NEB_5p--NEB_3p.bam; do 
	isoseq3 refine "$i" primers.fa "${i/.noPrimers_reads.NEB_5p--NEB_3p.bam/.FLNC_reads.bam}" --require-polya;
done
#Rename headers of fastq file with random SRR accessions. This step is necessary to run EGAPx
declare -A SAMPLES
SAMPLES=( ["Pupae.hifi.FLNC_reads.fastq"]="SRR9990001"
          ["Queens.hifi.FLNC_reads.fastq"]="SRR9990002"
          ["Workers.hifi.FLNC_reads.fastq"]="SRR9990003" )

for file in "${!SAMPLES[@]}"; do
    srr="${SAMPLES[$file]}"
    output="${file%.fastq}.srr.fastq"

    echo "Rinominando $file in $output con prefisso $srr..."

    awk -v srr="$srr" '
    BEGIN {count=0}
    {
        if (NR % 4 == 1) {
            count++;
            # Formato finale: @gnl|SRA|SRR999000X.count.1
            print "@gnl|SRA|" srr "." count ".1";
        }
        else if (NR % 4 == 3) {
            print "+";
        }
        else {
            print $0;
        }
    }' "$file" > "$output"
done
#Now we need to convert from fastq to fasta. Thesde files will be used for gene annotation
seqtk seq -A SRR9990001_reads.fastq > SRR9990001_reads.fasta
seqtk seq -A SRR9990002_reads.fastq > SRR9990002_reads.fasta
seqtk seq -A SRR9990003_reads.fastq > SRR9990003_reads.fasta
#Cluster reads (i.e isoforms)
for i in *FLNC_reads.bam; do 
	isoseq cluster2 --singletons --sort-threads 10 "$i" "${i/FLNC_reads.bam/FLNC.clusterd_reads.bam}";
done
