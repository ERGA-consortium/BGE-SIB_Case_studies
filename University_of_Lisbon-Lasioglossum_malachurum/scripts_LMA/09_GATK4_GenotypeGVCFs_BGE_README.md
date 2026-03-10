\# **Performs joint genotyping on the cohort-level GVCF using GATK4 GenotypeGVCFs** #





\\#### 1. Before Running ####



You must edit the script to set your own input/output directories and sample list path



\\\*Project\\\_path\\\* = \\\*<PATH/to/Project>\\\*



\\\*assembly\\\_path\\\* = \\\*<PATH/to/directory/masked\\\_reference\\\_genome>\\\*



\\\*assembly\\\* = <your\\\_reference\\\_genome\\\_masked>



\\\*cohort\\\_dir\\\* = \\\*<PATH/for/output\\\_cohort\\\_gvcf\\\_files>\\\*



---------------------------------------------------------------------------------



\\#### 2. Parameters used for bwa-mem ####



\\\*-R\\\* (reference genome)



\\\*-V\\\* (path to a list file containing the PATH to gvcf files)



\\\*-O\\\* (output directory)



---------------------------------------------------------------------------------



\\#### \\\*\\\*3. Activate the micromamba environment containing GATK4 GenotypeGVCFs\\\*\\\* ####



\\\*bash\\\*

\\\*micromamba activate <your\\\_environment\\\_name>\\\*



---------------------------------------------------------------------------------



\\#### 4. Run the script ####



\\\*bash 09\\\_GATK4\\\_GenotypeGVCFs\\\_BGE.sh\\\*



---------------------------------------------------------------------------------



