#!/bin/bash

# $1 is the assembly file to be evaluated

# Run QUAST to assess assembly quality
# -o sets the output directory name
# -t specifies the number of threads to use
quast.py -o quast_output_$1 $1 -t 8

