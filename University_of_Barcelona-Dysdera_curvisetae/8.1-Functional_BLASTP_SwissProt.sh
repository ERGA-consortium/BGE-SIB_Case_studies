#!/bin/bash

# $1 is the input protein FASTA file
input_fasta=$1

# Reference SwissProt protein database
database=/users-d1/silvia.garcia/Functional_annotation_final_version/SwissProt_03Aug2022/uniprot_sprot.fasta

echo "# Starting BLASTP..."
echo "# Input FASTA: ${input_fasta}"
date

# Run BLASTP against SwissProt
/users-d3/vadim.pisarenco/Programs_V/ncbi-blast-2.13.0+/bin/blastp \
  -query $input_fasta \
  -db $database \
  -out Dcurv_hap1_blastp.swissprot.outfmt6 \
  -evalue 1e-3 \
  -num_threads 48 \
  -outfmt "6 std qlen slen"

date
echo "# DONE"
