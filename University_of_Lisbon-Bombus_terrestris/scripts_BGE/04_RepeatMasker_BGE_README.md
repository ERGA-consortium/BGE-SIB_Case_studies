**# Runs RepeatMasker to soft-mask repetitive elements in a reference genome prior to alignment #**





\#### 1. Before Running ####



You must edit the script to set your own input/output directories



*species = Selects the appropriate repeat library*

*reference* = *reference genome used*



---------------------------------------------------------------------------------



\#### 2. Parameters used ####



*-pa 50* (number of threads used)

*-species* (library to be used)

*-xsmall* (Soft-mask repeats (lowercase))

*-gff* (Output annotation in GFF format)



---------------------------------------------------------------------------------



\#### **3. Activate the micromamba environment containing fastp** ####



*bash*

*micromamba activate <your\_environment\_name>*



---------------------------------------------------------------------------------



\#### 4. Run RepeatMasker ####



*bash* 

*RepeatMasker -pa 50 -species "<appropriate repeat library>" -xsmall -gff <reference>.fna*



---------------------------------------------------------------------------------









