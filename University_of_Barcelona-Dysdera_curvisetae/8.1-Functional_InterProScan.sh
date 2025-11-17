#!/bin/bash

# $1 is the input protein FASTA file
input_fasta=$1
output_fasta=$1.formatted.fasta

echo "# Starting InterProScan run..."
echo "# Input FASTA: ${input_fasta}"
echo "# Formatted FASTA (no '*' in sequences): ${output_fasta}"

# Remove '*' characters from sequences and save as new FASTA
sed 's/\*//g' $input_fasta > $output_fasta

date

# Run InterProScan
# -i        : input FASTA
# -t p      : input type proteins
# -iprlookup: include InterPro annotations
# -goterms  : include GO terms
# -dp       : disable pre-calculated match lookup (full run)
# -cpu 24   : number of threads
/users-d1/silvia.garcia/Genome_annotation/Functional_Annotation/Nesiotes_hap1/InterProScan_run/my_interproscan/interproscan-5.72-103.0/interproscan.sh \
  -i $output_fasta -t p -iprlookup -goterms -dp -cpu 24

date
echo "# DONE"

