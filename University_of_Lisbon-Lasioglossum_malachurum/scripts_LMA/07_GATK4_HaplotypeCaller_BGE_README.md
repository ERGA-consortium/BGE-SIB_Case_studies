\# **Runs GATK4 HaplotypeCaller per sample to produce GVCF files** #





\#### 1. Before Running ####



You must edit the script to set your own input/output directories and sample list path



\*Project\_path\* = \*<PATH/to/Project>\*



\*sample\_list\* = \*<PATH/to/sample\_list>\* (A text file containing one sample ID per line)



\*assembly\_path\* = \*<PATH/to/directory/masked\_reference\_genome>\*



\*assembly\* = <your\_reference\_genome\_masked>



\*bam\_dir\* = \*<PATH/to/output\_directory\_with\_dedup\_BAM\_files>\*



\*gvcf\_dir\* = \*<PATH/to/output\_directory\_for\_gvcf\_files>\*



---------------------------------------------------------------------------------



\#### 2. Parameters used for bwa-mem ####



\*-R\* (reference genome)



\*-I\* (deduplicated, sorted and indexed BAM files)



\*-ERC GVCF\* (enables GVCF mode. Per-sample gvcf file)



\*-O\* (output directory)



---------------------------------------------------------------------------------



\#### \*\*3. Activate the micromamba environment containing GATK4 HaplotypeCaller\*\* ####



\*bash\*

\*micromamba activate <your\_environment\_name>\*



---------------------------------------------------------------------------------



\#### 4. Run the script ####



\*bash 07\_GATK4\_HaplotypeCaller\_BGE.sh\*



---------------------------------------------------------------------------------



