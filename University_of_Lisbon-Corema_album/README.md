# Title of the case study
This repository includes all the scripts used to annotate the genome of *Corema album*

# Corema-album-annotation

## Bruno Nevado, 2025


# Notes: 
- Following (with a few changes) the annotation protocol here: https://begibmcp.notion.site/Genome-Annotation-Protocol-df7a2c4b05444290b879556efcd3b1e2, kindly provided by Aureliano Bombarely.  

- Annotation covered only the 13 largest scaffolds of the assembly, which represent the 13 chromosomes of *C. album*.  

- RNAseq data used from 16 samples, representing 2 individuals (1 male, 8 females), 4 tissues and 2 sampling times (early morning and late afternoon).  


## Repeat annotation

Software used:  
- repeatmodeler v 2.0.6  
- repeatmasker v 4.1.8  
- TEsorter v 1.5.0  

```
module load repeatmodeler/2.0.6

## 1- Generate the RepeatModeler sequence database
BuildDatabase --name coremaV1 Calbum.V1.fas

## 2- Run RepeatModeler2
RepeatModeler -database coremaV1 -threads 48 -LTRStruct

# Output:
#Program Time: 29:19:15 (hh:mm:ss) Elapsed Time
#Working directory:  /home/bnevado/corema/RM_3997133.MonMar171313582025
#may be deleted unless there were problems with the run.
#
#The results have been saved to:
#  /home/bnevado/corema/
#    coremaV1-families.fa  - Consensus sequences for each family identified.
#    coremaV1-families.stk - Seed alignments for each family identified.
#    coremaV1-rmod.log     - Execution log.  Useful for reproducing results.
#
#


## 3- Run RepeatMasker

module load repeatmasker/4.1.8
# Loading requirement: python/3.10.11 hmmer/3.2.1 trf/4.10.0
RepeatMasker -lib coremaV1-families.fa -xsmall -gff -e ncbi -pa 4 Calbum.V1.fas

## 4- Run TESorted

module load TEsorter/1.5.0
#  Loading requirement: python/3.10.11 hmmer/3.3.2 blast+/2.13.0
  
TEsorter coremaV1-families.fa -p 48 -db rexdb-plant
```

<Details>
<summary><span style="color:darkred">Repeat classification summary</span></summary>

- Repeat classification summary done with R package repeatR  
- In total, ca. 51% of the genome assembly consists of repeat elements, with the majority (30% of the assembly) being unknown.  

- Repeat annotation results with all classes.  
![Repeat content per class, including unknown elements](./annotationResults/repearContentPerClass.1.png){#id .class width=60% height=60%}
- Repeat annotation results excluding unclassified repeats.  
![Repeat content per class, excluding unknown elements](./annotationResults/repearContentPerClass.2.png){#id .class width=60% height=60%}

</Details>

## Evidence-based annotation with RNAseq data

Software used:  
- trimgalore v 0.6.10  
- fastqc v 0.12.1  
- STAR v 2.7.11b  
- samtools/bcftools v1.17  
- stringtie v 3.0.1  
- transdecoder v 5.7.1
- blast+ v 2.13.0  
- diamond v 2.1.11  
- hmmer v 3.3.2
```
# 1. Read trimming and adaptor removal
module load trimgalore/0.6.10
module load fastqc/0.12.1
