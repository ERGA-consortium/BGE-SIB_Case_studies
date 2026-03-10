#! /bin/bash

#### Script to run fastqc e multiqc nos raw_data  ### amostras BGE

#### Directoria Operacional do Projecto   ##data3
Project_dir=/mnt/data3/Project_BGE_SSF_LM

####  Directoria da lista de amostras do projecto
samplelist=$Project_dir/Operational/2_metadata/sample_list_Lmal.txt
rawfastqsuffix1=_1.fastq.gz
rawfastqsuffix2=_2.fastq.gz


#### Directoria de input dos dados
inputdir=/Operational/4_data_analysis/2_fastp/

#### Directoria de output dos resultados fastqc
outputdir=/Operational/4_data_analysis/1_qc/2_fastqc_posfastp/


#### Directoria de output dos resultados multiqc
outputdir_multiqc=/Operational/4_data_analysis/1_qc/2_fastqc_posfastp/




for sample in `cat  $samplelist`; do


    fastqc -t 8  $Project_dir$inputdir$sample$rawfastqsuffix1  $Project_dir$inputdir$sample$rawfastqsuffix2  -o $Project_dir$outputdir

     
done


    multiqc  $Project_dir$outputdir  --outdir  $Project_dir$outputdir_multiqc  

