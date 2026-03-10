**# Runs fastp on raw paired-end data for trimming and filtering #**





\#### 1. Before Running ####



You must edit the script to set your own paths



*Project\_dir* = <*PATH/to/Project\_directory>*

*samplelist* = <*PATH/to/sample\_list>* (A text file containing one sample ID per line)

*rawfastqsuffix1* = *<sample>\_1.fastq.gz* (Expected file suffix for Read\_1)

*rawfastqsuffix2* = *<sample>\_2.fastq.gz* (Expected file suffix for Read\_2)

*inputdir* = <*PATH/to/raw\_data>* (A folder containing raw data)

*outputdir\_fastp* = <*PATH/to/outpudir\_for\_fastp>*



---------------------------------------------------------------------------------



\#### 2. Parameters used ####



*-g* (Trim Poly G)

*-l 50* (reads shorter than 50bp are removed)

*-q 20* (the quality value that a base is qualified. Value of 20 means phred quality >=Q20 is qualified)

*-u 30* (discard read if >30 % bases fall below Q20)

*--html* (generates .HTML reports)

*-j* (generates .JSON reports)

*--failed\_out* (directory to store the reads that failed)

---------------------------------------------------------------------------------



\#### **3. Activate the micromamba environment containing fastp** ####



*bash*

*micromamba activate <your\_environment\_name>*



---------------------------------------------------------------------------------



\#### 4. Run the script ####



*bash 02\_fastp\_rawdata\_BGE.sh*



---------------------------------------------------------------------------------











