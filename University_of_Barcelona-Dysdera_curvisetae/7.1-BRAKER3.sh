#!/bin/bash

# $1 is the assembly (masked genome) to annotate
# $2 is the working directory for BRAKER output

ASSEMBLY=$1
wd=$2

# Run BRAKER3 for genome annotation using RNA-seq and protein evidence
# All output and errors are redirected to braker3_output.log

time braker.pl \
  --genome=${ASSEMBLY}.fasta.masked \
  --prot_seq=/users-d1/silvia.garcia/Genome_annotation/DsilJoelProt_ArthODB10.fasta \
  --rnaseq_sets_ids=HHNM5DSXF_1_334UDI-idt-UMI_paired,HHNM5DSXF_1_346UDI-idt-UMI_paired,HHNM5DSXF_2_334UDI-idt-UMI_paired,HHNM5DSXF_2_346UDI-idt-UMI_paired,HHNM5DSXF_3_334UDI-idt-UMI_paired,HHNM5DSXF_3_346UDI-idt-UMI_paired,HHNM5DSXF_4_334UDI-idt-UMI_paired,HHNM5DSXF_4_346UDI-idt-UMI_paired \
  --rnaseq_sets_dirs=/users-d1/silvia.garcia/cnag_2024_2025/trimmomatic_results/trimmomatic_results_15_paired/ \
  --gff3 \
  --workingdir=$wd \
  --threads 48 \
  &> braker3_output.log
