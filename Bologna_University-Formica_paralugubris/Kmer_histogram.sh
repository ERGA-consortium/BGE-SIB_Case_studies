#!/bin/bash

READS=$1	#Reads
PREFIX=$2	#Prefix of output files
mkdir -p ../../Analyses/KMER_Hist

#Produce kmer histograms for different k-mer sizes
declare -a arr=("19" "21" "23" "25" "27" "29" "31")
for i in "${arr[@]}"
do
   kat hist --mer_len "$i"  --high 5000000 --hash_size 1000000000 --threads 10 -o ../../Analyses/KMER_Hist/Fpar_"$PREFIX"_Kmer"$i".hist "$READS"
done
