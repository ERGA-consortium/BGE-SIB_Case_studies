\# **Applies hard filters to raw cohort VCF using GATK4 VariantFiltration** #





\\#### 1. Before Running ####



You must edit the script to set your own input/output directories



\\\*Project\\\_path\\\* = \\\*<PATH/to/Project>\\\*



\\\*assembly\\\_path\\\* = \\\*<PATH/to/directory/masked\\\_reference\\\_genome>\\\*



\\\*assembly\\\* = <your\\\_reference\\\_genome\\\_masked>



\\\*vcf\\\_dir\\\* = \\\*<PATH/to/directory\\\_with\\\_cohort\\\_gvcf\\\_files>\\\*



---------------------------------------------------------------------------------



\\#### 2. Parameters used for bwa-mem ####



\\\*-R\\\* (reference genome)



\\\*-V\\\* (input raw cohort vcf)



\\\*--filter-expression "QD < 2.0 || FS > 60.0 || MQ < 40.0"\\\* (filters based on quality filters)



\\\*--filter-name "basic\_snp\_filter"\\\* (name assigned to filtered variants)



\\\*-O\\\* (output directory)



---------------------------------------------------------------------------------



\\#### \\\*\\\*3. Activate the micromamba environment containing GATK4 VariantFiltration\\\*\\\* ####



\\\*bash\\\*

\\\*micromamba activate <your\\\_environment\\\_name>\\\*



---------------------------------------------------------------------------------



\\#### 4. Run the script ####



\\\*bash 10\\\_GATK4\\\_VariantFiltration\\\_BGE.sh\\\*



---------------------------------------------------------------------------------



