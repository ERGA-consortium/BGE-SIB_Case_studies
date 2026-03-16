#!/bin/bash

GENOME=$1
READS=$2

#Produce kmer spectra
kat comp  -o ../../Analyses/Assembly/Default/Fpar_KmerSpectra -t 10 "$READS" "$GENOME"

#Busco
busco -i "$GENOME" -l insecta -o ../../Analyses/Assembly/Default/Fpar_Inse -m geno -c 5
busco -i "$GENOME" -l hymenoptera -o ../../Analyses/Assembly/Default/Fpar_Inse -m geno -c 5
