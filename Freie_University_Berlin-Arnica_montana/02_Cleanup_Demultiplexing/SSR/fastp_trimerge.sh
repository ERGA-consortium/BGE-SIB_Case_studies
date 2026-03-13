#!/bin/bash

echo "Batch trim & merge Illumina read pairs with fastp."
echo "Usage: fastp_trimerge Inputdir Outputdir Suffix"
echo "Outputdir: path relative to directory above Input"

echo "Input dir: $1"
echo "Output dir: $2"
echo "Suffix: $3"
curdir=$(pwd)

cd $1
mkdir ../$2
mkdir ../${2}/report

for fwread in *_R1${3}; do
	rvread=${fwread%_R1${3}}_R2${3}
	outfile=${fwread%_R1${3}}_M${3}
	fastp --in1 ${fwread} --in2 ${rvread} --merge --merged_out ../$2/${outfile} --detect_adapter_for_pe --html ../$2/report/${fwread%_R1${3}}_fastp.html --json ../$2/report/${fwread%_R1${3}}_fastp.json
done
cd $curdir