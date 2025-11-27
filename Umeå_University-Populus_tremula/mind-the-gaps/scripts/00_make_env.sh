#!/bin/bash -l 

# make the environment
mamba create -n chromcomp python=3.11 -y

# activate 
mamba activate chromcomp

# install packages
mamba install -c bioconda \
	bwa=0.7.19 \
	samtools=1.21 \
	gatk4=4.6.2.0 \
	plink2=2.00a3 \
	vcftools=0.1.17 \
	bcftools=1.21 \
	busco=5.8.2 \
	repeatmasker=4.1.8 \
	fastqc=0.12.1 \
	multiqc=1.29 \
	fastp=1.0.1 \
	-y

mamba install -c conda-forge \
	openjdk>=17    
	

mamba install -c bioconda \
    nextflow=25.04.2
