**# Runs FastQC/MultiQC on filtered paired-end data #**





\#### 1. Before Running ####



You must edit the script to set your own input/output directories and sample list path



*Project\_dir* = *<PATH/to/Project\_directory>*

*samplelist* = *<PATH/to/sample\_list>* (A text file containing one sample ID per line)

*rawfastqsuffix1* = *<sample>\_1.fastq.gz* (Expected file suffix for Read\_1)

*rawfastqsuffix2* = *<sample>\_2.fastq.gz* (Expected file suffix for Read\_2)

*inputdir* = *<PATH/to/raw\_data>* (A folder containing raw data)

*outputdir\_fastQC* = *<PATH/to/outpudir\_for\_fastQC>*

*outputdir\_multiQC* = *<PATH/to/outpudir\_for\_multiQC>*



---------------------------------------------------------------------------------



\#### 2. Activate the micromamba environment containing FastQC and MultiQC ####



*bash*

*micromamba activate <your\_environment\_name>*



---------------------------------------------------------------------------------



\#### 3. Run the script ####



*bash 03\_fastqc\_multiqc\_filtered\_data\_BGE.sh*



---------------------------------------------------------------------------------

