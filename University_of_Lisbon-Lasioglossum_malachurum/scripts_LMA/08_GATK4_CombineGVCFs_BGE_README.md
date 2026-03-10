\# **Combines multiple per-sample gvcfs into a single cohort-level gvcf using GATK4 CombineGVCFs** #



\\#### 1. Before Running ####







You must edit the script to set your own input/output directories



\\\*Project\\\_path\\\* = \\\*<PATH/to/Project>\\\*



\\\*gvcf\\\_dir\\\* = \\\*<PATH/to/directory\\\_with\\\_gvcf\\\_files>\\\*



\\\*output\\\_dir\\\* = \\\*<PATH/to/output\\\_directory>\\\*



\\\*assembly\\\_path\\\* = \\\*<PATH/to/directory/masked\\\_reference\\\_genome>\\\*



\\\*assembly\\\* = <your\\\_reference\\\_genome\\\_masked>





---------------------------------------------------------------------------------



\\#### 2. Parameters used for bwa-mem ####



\\\*-R\\\* (reference genome)



\\\*-V\\\* (path to a list file containing the PATH to gvcf files)



\\\*-O\\\* (output directory)



---------------------------------------------------------------------------------



\\#### \\\*\\\*3. Activate the micromamba environment containing GATK4 CombineGVCFs\\\*\\\* ####



\\\*bash\\\*

\\\*micromamba activate <your\\\_environment\\\_name>\\\*



---------------------------------------------------------------------------------



\\#### 4. Run the script ####



\\\*bash 08\\\_GATK4\\\_CombineGVCFs\\\_BGE.sh\\\*



---------------------------------------------------------------------------------



