# **Align fastp-filtered reads using bwa-mem, followed by sorting, indexing, and QC statistics using samtools** #


#### 1. Before Running ####

You must edit the script to set your own input/output directories and sample list path

*Project_path* = *<PATH/to/Project>*

*sample_list* = *<PATH/to/sample_list>* (A text file containing one sample ID per line)

*fastq_path* = *<PATH/to/filtered_fastp_files>*

*assembly_path* = *<PATH/to/directory/masked_reference_genome>*

*output_base* = *<PATH/to/output_directory>*

*bam_dir* = *<PATH/to/output_directory_for_BAM_files>*

*qc_dir* = *<PATH/to/output_directory_for_qc_files>*

*assembly* = <your_reference_genome_masked>

---------------------------------------------------------------------------------

#### 2. Parameters used for bwa-mem ####

*-M* (mark shorter slipt hits as secondary)

*-t* (number of threads used)

*-R* (Read group information; adds metadata to every read)

*-gff* (Output annotation in GFF format)

---------------------------------------------------------------------------------

#### **3. Activate the micromamba environment containing bwa-mem+samtools** ####

*bash*
*micromamba activate <your_environment_name>*

---------------------------------------------------------------------------------

#### 4. Run the script ####

*bash 05_bwamem_samtools_BGE.sh*



---------------------------------------------------------------------------------



























