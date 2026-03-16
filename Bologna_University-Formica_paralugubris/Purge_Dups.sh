#!/bin/bash

READS=$(realpath $1)
GENOME=$(realpath $2)

mkdir ../../Analyses/Assembly/Default/Purge_dups
cd ../../Analyses/Assembly/Default/Purge_dups

minimap2 -xasm20 -t 10 "$GENOME" "$READS" | gzip -c - > Fpar.paf.gz
pbcstat Fpar.paf.gz
calcuts PB.stat > cutoffs 2> calcults.log

split_fa "$GENOME" > Fpar.split
minimap2 -xasm5 -DP Fpar.split Fpar.split | gzip -c - > Fpar.split.self.paf.gz
purge_dups -2 -T cutoffs -c PB.base.cov Fpar.split.self.paf.gz > dups.bed 2> purge_dups.log
get_seqs -e dups.bed "$GENOME"
