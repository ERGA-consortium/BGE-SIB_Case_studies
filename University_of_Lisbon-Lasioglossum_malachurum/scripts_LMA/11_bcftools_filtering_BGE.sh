#! /bin/bash

#### Script to run bcftools on vcf files output of freebayes  ### bcftools filter

#### vcf file de input
vcffile=cohort.GCA_028455605.1_LMAL_v2.2.0_genomicCHROM.masked.flt.vcf.gz

#### genoma de referecia
ref_genome=/mnt/data3/Project_BGE_SSF_LM/Operational/3_reference_genomes/RepeatMasker_CHROM/GCA_028455605.1_LMAL_v2.2.0_genomicCHROM.masked.fna

snp=snps_
Q20=Q20_
norm=norm_
bi=bi_
GQ=GQ20_
emiss08=emiss08_  ##excluir mais de 0.8 missing data
emiss05=emiss05_  ##excluir mais de 0.5 missing data
mis=mis_

### file de output de contagens
outtxt=out_samples_BGE_LM_v1.txt



## variação total

zgrep -v -c '^#' $vcffile > $outtxt

##snps
bcftools view $vcffile -i 'TYPE="snp"' -Oz -o $snp$vcffile

zgrep -v -c '^#' $snp$vcffile >> $outtxt


##norm
bcftools norm --fasta-ref $ref_genome -m -snps $snp$vcffile -Oz -o $snp$norm$vcffile

zgrep -v -c '^#' $snp$norm$vcffile >> $outtxt

## biallelicos
bcftools view $snp$norm$vcffile --min-alleles 2 --max-alleles 2  -Oz -o  $snp$norm$bi$vcffile  

zgrep -v -c '^#' $snp$norm$bi$vcffile >> $outtxt


##Qual>20
bcftools filter -i 'QUAL>=20' $snp$norm$bi$vcffile  -Oz -o  $snp$norm$bi$Q20$vcffile

zgrep -v -c '^#' $snp$norm$bi$Q20$vcffile >> $outtxt

##GQ<20

bcftools filter  --exclude 'GQ<20' --set-GTs '.' $snp$norm$bi$Q20$vcffile -Oz -o $snp$norm$bi$Q20$GQ$vcffile

zgrep -v -c '^#' $snp$norm$bi$Q20$GQ$vcffile >> $outtxt


## exluir misisng 0.8
bcftools view  $snp$norm$bi$Q20$GQ$vcffile -e 'F_MISSING>0.8'  -Oz -o $snp$norm$bi$Q20$GQ$emiss08$vcffile

zgrep -v -c '^#' $snp$norm$bi$Q20$GQ$emiss08$vcffile >> $outtxt

## exluir misisng 0.5
bcftools view  $snp$norm$bi$Q20$GQ$vcffile -e 'F_MISSING>0.5'  -Oz -o $snp$norm$bi$Q20$GQ$emiss05$vcffile

zgrep -v -c '^#' $snp$norm$bi$Q20$GQ$emiss05$vcffile >> $outtxt


## exluir missing
bcftools view -e 'GT[*] = "mis"'  $snp$norm$bi$Q20$GQ$vcffile -Oz  -o $snp$norm$bi$Q20$GQ$mis$vcffile

zgrep -v -c '^#' $snp$norm$bi$Q20$GQ$mis$vcffile >> $outtxt