# Raw Data Processing

3RAD data are quadruple indexed, i.e. the combination of outer indices and inner barcodes identifies each sample. In our case, only the P7 index identifies libraries (and should thus be used for an initial step of demultiplexing), as the P5 index is used to distinguish PCR duplicates. Thus, the indices need to be read out and stored as well as the actual insert / "read" part of the library constructs.

## Prerequisites
Beside the file path to the actual data (run folder, name starts with date; hard coded in script), this bash script needs the ["bcl2fastq" software](https://support.illumina.com/content/dam/illumina-support/documents/documentation/software_documentation/bcl2fastq/bcl2fastq_letterbooklet_15038058brpmi.pdf) provided by Illumina, and is supposed to be run on an HPC with the SLURM task manager.
