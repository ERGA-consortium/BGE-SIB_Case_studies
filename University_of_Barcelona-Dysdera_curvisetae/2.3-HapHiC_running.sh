#!/bin/bash

# $1 is the assembly file used as the reference genome

# Run the HapHiC pipeline to phase and scaffold the assembly using Hi-C data
# Arguments:
#   ../$1                     : path to the reference assembly (relative path)
#   ../HiC_$1_filtered.bam    : filtered Hi-C alignment file
#   11                        : number of threads or processing cores to use
#   --RE "GATC,TTAA"          : restriction enzyme recognition sites used for Hi-C library preparation
#   --correct_nrounds 2       : number of rounds of contact correction to improve phasing accuracy
/users-d1/silvia.garcia/Programs/HapHiC/haphic pipeline ../$1 ../HiC_$1_filtered.bam 11 --RE "GATC,TTAA" --correct_nrounds 2

