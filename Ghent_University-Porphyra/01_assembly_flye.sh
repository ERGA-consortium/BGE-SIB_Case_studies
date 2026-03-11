#!/bin/bash
# Assembly for ONT reads
# Param: $1 = raw_ont_reads.fastq.gz
flye --nano-raw $1 --genome-size 100m --out-dir ./01_flye_output --threads 16
