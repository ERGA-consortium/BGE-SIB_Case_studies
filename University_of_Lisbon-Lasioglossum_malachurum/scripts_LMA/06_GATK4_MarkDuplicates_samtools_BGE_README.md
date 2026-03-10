\# **Marks PCR duplicates in BAM files using GATK4 MarkDuplicates** #

\# **Indexes the resulting deduplicated BAM files using samtools index** #





\#### 1. Before Running ####



You must edit the script to set your own input/output directories and sample list path



\*Project\_path\* = \*<PATH/to/Project>\*



\*sample\_list\* = \*<PATH/to/sample\_list>\* (A text file containing one sample ID per line)



\*assembly\_path\* = \*<PATH/to/directory/masked\_reference\_genome>\*



\*assembly\* = <your\_reference\_genome\_masked>



\*output\_base\* = \*<PATH/to/output\_directory>\*



\*bam\_dir\* = \*<PATH/to/output\_directory\_for\_BAM\_files>\*



\*qc\_dir\* = \*<PATH/to/output\_directory\_for\_qc\_files>\*



---------------------------------------------------------------------------------



\#### 2. Parameters used for bwa-mem ####



\*-I\* (input BAM files)



\*-O\* (output dedup BAM files)



\*-M\* (generates a metrics file)



\*--REMOVE\_DUPLICATES true/false\* (GATK recommends to mark them, not to remove them)



---------------------------------------------------------------------------------



\#### \*\*3. Activate the micromamba environment containing GATK4 MarkDuplicates +samtools\*\* ####



\*bash\*

\*micromamba activate <your\_environment\_name>\*



---------------------------------------------------------------------------------



\#### 4. Run the script ####



\*bash 06\_GATK4\_MarkDuplicates\_samtools\_BGE.sh\*



---------------------------------------------------------------------------------



