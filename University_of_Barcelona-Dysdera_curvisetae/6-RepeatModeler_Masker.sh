#!/bin/bash

# $1 is the input FASTA file (without extension)

FASTAFILE=$1

# Build a repeat library database
/users-d1/silvia.garcia/Repetitive_Elements/RepeatModeler/RepeatModeler-2.0.3/BuildDatabase -name $FASTAFILE $FASTAFILE.fasta

# Run RepeatModeler to find repeats
/users-d1/silvia.garcia/Repetitive_Elements/RepeatModeler/RepeatModeler-2.0.3/RepeatModeler -database $FASTAFILE -pa 36 -LTRStruct

# Get the latest output directory from RepeatModeler
LAST_DIR=$(ls -dt */ | head -n 1)

# Mask repeats in the assembly using RepeatMasker with the custom library
/users-d1/silvia.garcia/Repetitive_Elements/RepeatMasker/RepeatMasker/RepeatMasker $FASTAFILE.fasta -pa 36 -gff -xsmall -s -gccalc -a -engine ncbi -lib $LAST_DIR/consensi.fa.classified -dir MaskerOutput_$FASTAFILE

