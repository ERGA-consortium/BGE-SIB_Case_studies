#!/bin/bash

READS=$1

mkdir ../../Analyses/Assembly/Deafult
mkdir ../../Analyses/Assembly/Agressive_Purging
mkdir ../../Analyses/Assembly/Agressive_Purging_2

#Default parameters
hifiasm -o ../../Analyses/Assembly/Deafult/Fpar_HiFiasm_NoCont -t15 "$READS"
#More agressive purging
hifiasm -o ../../Analyses/Assembly/Agressive_Purging/Fpar_HiFiasm_NoCont_Agg.Purg -s 0.45 -t 15 "$READS"
#More agressive purging 2
hifiasm -o ../../Analyses/Assembly/Agressive_Purging_2/Fpar_HiFiasm_NoCont_Agg.Purg -s 0.35 -t 15 "$READS"
