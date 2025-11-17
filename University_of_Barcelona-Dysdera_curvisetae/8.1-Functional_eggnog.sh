#!/bin/bash

# $input_fasta  : input protein FASTA file
# $output       : base name for output files
# $ncpu         : number of CPU threads

echo "# Starting eggNOG run..."
echo "# Input FASTA: ${input_fasta}"
date

# 1) Search step: find orthologs using DIAMOND, without annotation
emapper.py -m diamond --sensmode ultra-sensitive --no_annot \
            -i $input_fasta --itype proteins \
            -o $output --cpu $ncpu

# 2) Annotation step: annotate the previously found orthologs
emapper.py -m no_search --annotate_hits_table $output.emapper.seed_orthologs \
            -o "${output}_annot_1" --dbmem --cpu $ncpu

# 3) Refined annotation targeting specific taxa and reporting one-to-one orthologs
# Taxids:
#   Arthropoda  6656
#   Chelicerata 6843
#   Arachnida   6854
emapper.py -m no_search --annotate_hits_table $output.emapper.seed_orthologs \
            -o "${output}_annot_2" --dbmem --cpu $ncpu \
            --report_orthologs --target_orthologs one2one \
            --target_taxa 6854,6843,6656

date
echo "# DONE"
