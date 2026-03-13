# SSRseq cleanup
As the SSRseq data already comes demultiplexed by sample (based on the p5 and p7 indexes) and will be filtered based on the contained F+R primer sequences during genotyping, this step simply consists of quality filtering and read pair merging with fastq.

## Prerequisites
Beside the demultiplexed raw data, this script needs an installation of [fastp](https://github.com/OpenGene/fastp) and can be run on a laptop / desktop computer. The usage of the script is:
  
  fastp_trimerge Inputdir Outputdir Suffix

Where Outputdir is a path relative to the directory above the Inputdir, and Suffix is the identical file ending after "R1" (forward reads) or "R2" (reverse reads).
