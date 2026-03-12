# Hybridization and its implication for mountain hares
## BGE Project 11; University of Eastern Finland


Script collection for the L. timidus and L. europaeus population genomics BGE project

**From alignment to analysis**

**Warning!** This is not a dev project. This is a documentation of a tailored analysis pipeline, for our cluster and our goals.
For rerunning, you might need to adjust parameters for your computer/cluster.

Note: LE and LT versions of the same script are usually the same especially if/when not both are included here. They were just used to run both species' processing in parallel.

### Data

ENA accession for (demultiplexed) fastq.gz files: PRJEB101756


### Reference genomes

Related reference genomes were assembled under different projects  
Their accessions and connected publications are available here:  

| Species	       | NCBI GenBank    | DOI                    |
-----------------------|-----------------|------------------------|
| _Lepus europaeus_    | GCA_033115175.1 | 10.24072/pcjournal.393 |
| _Lepus timidus_      | GCA_033115175.1 | 10.24072/pcjournal.514 |


For the purpose of this project, we limited analysis to using the _L. europaeus_ genome as a reference.


### Environment

All processing and analysis was run on the CSC cluster of Finland. Script parameters are adjusted to this.  
Software versions are included in the scripts (specified at loading the relevant modules)

### Metadata

Metadata files include:
- sample lists
- sampling information (SAMPLE_metadata.txt)
- sex identification based on allosome read coverage (Hybrid.sex_and_Y.v0903.tsv)

### Quality control and demultiplexing

- *fqc.sh* : FastQC run for all samples
- *demult.LE.sh* : Demultiplexing _L. europaeus_ data
- *demult.LT.sh* : Demultiplexing _L. timidus_ data

FastQC was run both before and after demultiplexing


### Alignment

- *bwa_index.sh* : indexing the reference genome
- *bwa_align.sh* : aligning samples to reference genome

Sample lists are provided in the Metadata/*.list files


### Variant calling and filtering

#### Variant calling

- *hybrid_varcall.autosomes.LE.sh* : Autosome (diploid) variant calling
- *radseq_femaleX_varcall.sh* : Female X chromosome (diploid) variant calling
- *radseq_haploid_varcall.sh* : Haploid variant calling for male X and Y chromosomes


#### Genotyping and quality filtering

- *genotype.all.array.sh* : Combining and joint genotyping called variants, paralellized for chromosomes
- *vcf_combine_Y.sh* : While variant calling and combining the Y g.vcf files was done separately to account for haploidy, further filters were the same
- *radseq_vcf_process.sh* : Genotype quality filtering across all samples combined

### Analysis

Note, that most analysis steps were run in interactive sessions, so that we could easily adjust parameters and rerun where needed.
Any code provided here is for documenting analysis parameters

Analysis steps include:
- PCA (plink v2.00a5)
- Interclass heterozygosity with triangulaR (triangulaR v0.0.1) and vcfR (vcfR v1.15.0)
- ADMIXTURE (ADMIXTURE v1.3.0)
- Fst and Tajima-D calculation (vcftools 0.1.17)
- ELAI run (ELAI v1.01) + summary and mean introgression values calculated in R (base R v4.4.3, tidyverse v2.0.0)

R scripts are mainly for plotting, with light postprocessing included


#Authors

Jaakko Pohjoismäki (https://orcid.org/0000-0002-1185-3610), University of Eastern Finland, Finland 
Zsófia Fekete (https://orcid.org/0000-0002-9086-5459), University of Oulu, Oulu, Finland


# Acknowledgements

This work was supported by the SIB Swiss Institute of Bioinformatics under the Biodiversity Genomics Europe project.