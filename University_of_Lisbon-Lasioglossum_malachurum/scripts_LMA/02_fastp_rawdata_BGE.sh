#! /bin/bash

#### Script to run fastp nos rawdata  ### amostras  BGE


#### Directoria Operacional do Projecto   ##data3
Project_dir=/mnt/data3/Project_BGE_SSF_LM

####  Directoria da lista de amostras do projecto
samplelist=$Project_dir/Operational/2_metadata/sample_list_Lmal.txt
rawfastqsuffix1=_1.fastq.gz
rawfastqsuffix2=_2.fastq.gz


#### Directoria de input dos dados  ###data1 
inputdir=/mnt/data1/Project_BGE_SSF_LM/Operational/1_raw_data/

#### Directoria de output dos resultados fastp
outputdir=/Operational/4_data_analysis/2_fastp/



for sample in `cat  $samplelist`; do


       fastp \
       -i $inputdir$sample$rawfastqsuffix1 \
       -I $inputdir$sample$rawfastqsuffix2 \
       -o $Project_dir$outputdir$sample$rawfastqsuffix1 \
       -O $Project_dir$outputdir$sample$rawfastqsuffix2 \
       -g -l 50 -q 20 -u 30  \
       --html $Project_dir$outputdir$sample'.html' \
       -j $Project_dir$outputdir$sample'.json' \
       --failed_out $Project_dir$outputdir$sample'.txt'
                
done


### -u discard read if >30 % bases fall below Q20
### -g -l 80 -p -D 5 -q 20
### -g trim ply G
### -l 80   read shorter than are removed
### -p   -overrepresentation_analysis    enable overrepresented sequence analysis. 
### -q  --qualified_quality_phred      the quality value that a base is qualified. Default 15 means phred quality >=Q15 is qualified. (int [=15])
### -D 5 depup