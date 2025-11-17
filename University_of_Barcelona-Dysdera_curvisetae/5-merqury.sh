#!/bin/bash

# $1 is the assembly file

# Variables
ASSEMBLY=$1                            # Assembly to evaluate
READS_HIFI="JR3_CUR_35.hifi_reads.filt.fastq.gz"  # PacBio HiFi reads
K=25                                   # k-mer size (recommended: 21 or 23)
THREADS=24                             # Number of CPU threads
OUTDIR=$1"_merqury_out_k${K}"         # Output directory

# Count k-mers in the HiFi reads using Meryl
meryl count k=${K} threads=${THREADS} output $1_reads_hifi.k${K}.meryl ${READS_HIFI}

# Run Merqury to assess assembly quality using the k-mer database
merqury.sh $1_reads_hifi.k${K}.meryl ${ASSEMBLY} ${OUTDIR}
