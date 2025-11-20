#!/bin/bash

# $1 is the sample name / assembly identifier

##
# Additional step: TSEBRA to refine gene models
##
# BRAKER3 can be very strict in some cases, potentially missing real genes.
# Following the developers' recommendations, we use TSEBRA to combine evidence 
# from Augustus and GeneMark-ETP and keep supported models.

# Run TSEBRA combining Augustus hints and GeneMark-ETP supported predictions
tsebra.py -g /users-d1/silvia.garcia/Genome_annotation_final_version/$1/braker3_output/Augustus/augustus.hints.gtf,/users-d1/silvia.garcia/Genome_annotation_final_version/$1/braker3_output/GeneMark-ETP/genemark_supported.gtf \
          -k /users-d1/silvia.garcia/Genome_annotation_final_version/$1/braker3_output/GeneMark-ETP/genemark_supported.gtf \
          -e /users-d1/silvia.garcia/Genome_annotation_final_version/$1/braker3_output/hintsfile.gff \
          -o tsebra_both_keep_genemark.gtf 2> tsebra.log

# Extract the longest isoform per gene
get_longest_isoform.py -g tsebra_both_keep_genemark.gtf -o tsebra_both_keep_genemark_longest.gtf

# Generate FASTA files for transcripts and proteins from the refined annotation
getAnnoFastaFromJoingenes.py -g $1 \
                              -o tsebra_both_keep_genemark_longest \
                              -f tsebra_both_keep_genemark_longest.gtf

