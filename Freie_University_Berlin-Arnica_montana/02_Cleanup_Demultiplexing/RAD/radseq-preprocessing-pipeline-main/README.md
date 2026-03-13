# RADseq preprocessing pipeline

User friendly automated pipeline for preprocessing RAD Seq data using Snakemake.

The following is an instruction on how to run the preprocessing pipeline for *ddRAD-like* data.  The pipeline will automatically execute the following steps:

**1. PhiX filter**

* Goal: Remove possible contamination with Enterobacteria phage phiX174 (in case you added phiX control library to your sequencing).
* Strategy: pooled library paired-end reads are mapped to the Enterobacteria phage phiX174 reference genome (NC_001422.1) using [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml) [1] with default parameters, and the unmapped reads are saved for the following processing steps.

**2. Demultiplexing**

* Goal: demultiplex reads based on inline barcodes and assign them to their correct sample ID
* Strategy: use the software [Flexbar](https://github.com/seqan/flexbar/wiki/Manual) [2]. After demultiplexing, the pipeline renames the individuals.

**3. PCR duplicates filter**

* Goal: remove PCR duplicates based on recognition of identical reads with the same iTru-8N index sequence.
* Strategy: use a custom python script that outputs only one copy of each unique DNA molecule. After this step, the pipeline concatenated the replicates of each library.

**4. Merge paired-end reads together and plot read lengths of these merged reads**

* Goal: merge forward and reverse reads using [PEAR](https://cme.h-its.org/exelixis/web/software/pear/doc.html) v.0.9.11 [3].
* Strategy: use software Pear to merge R1 and R2. The Snakemake pipeline generates a series of plots with the fragment length distribution of each individual.

**5. *In Silico* digestion of paired reads**

* Goal: recognize and filter out undigested and chimeric sequences.
* Strategy: use a custom python script

**6. Restriction enzyme site filtering**

* Goal: filter out reads that don't start and end with the correct restrictions sites.
* Strategy: use a custom python script that outputs only reads with correct sequences at both ends.

**7. Quality score filtering**

* Goal: filter out low quality (phred score lower than 30) reads.
* Strategy: use the software [Trimmomatic](http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/TrimmomaticManual_V0.32.pdf) [4]

:memo:**References**

- [1] Langmead, B., & Salzberg, S. L. (2012). Fast gapped-read alignment with Bowtie 2. Nature Methods, 9(4), 357–359.
- [2] Roehr, J. T., Dieterich, C., & Reinert, K. (2017). Flexbar 3.0 – SIMD and multicore parallelization. Bioinformatics , 33(18), 2941–2942.
- [3] Zhang, J., Kobert, K., Flouri, T., & Stamatakis, A. (2014). PEAR: a fast and accurate Illumina Paired-End reAd mergeR. Bioinformatics , 30(5), 614–620.
- [4] Bolger, A. M., Lohse, M., & Usadel, B. (2014). Trimmomatic: a flexible trimmer for Illumina sequence data. Bioinformatics , 30(15), 2114–2120.


With just a few simple alterations to the configuration files, the user can take advantage of the modular, portable, and scalable capabilities of Snakemake - including but not limited to the very useful automatic software dependency resolution.  

This workflow is designed to work on Linux-based operating systems.  HPC Cluster execution is not currently possible.

---

## Instructions

1. Downloading the workflow folder
1. Installing Conda
1. Creating our Snakemake radPreprocessing environment
1. Modifying our configuration files
1. Running the workflow

---

### Downloading the workflow


To clone the repository, use the following command:

```
git clone https://git.imp.fu-berlin.de/begendiv/radseq-preprocessing-pipeline.git
```

---

### Installing Conda

- Conda (v4.11+)  *but may work on older versions*

If you already have conda installed on your system, please skip to [Creating our Snakemake conda environment](#creating-our-snakemake-conda-environment)

Download the linux Miniconda3 installer from the following URL: https://docs.conda.io/en/latest/miniconda.html

Run the miniconda3 installation and check if it worked:

```
bash /<your_path_to>/Miniconda3-latest-Linux-x86_64.sh
##Follow miniconda3 installation instructions##

source ~/.bashrc

conda update conda
```

If  `conda command not found` please close and re-open your terminal for conda installation to take effect, and then update.

---
<div id="creating-our-snakemake-conda-environment"></div>

## Creating our RADseq preprocessing (Snakemake) Conda environment

The pipeline requires the following software to run:

- snakemake = `7.25.0`
- python >= `3.10`
- tabulate = `0.8.10`
- mamba = `1.4.1`
- pandas = `1.4.2`

The easiest method to install this software stack is to create a RADseq preprocessing conda environment with the provided `snakemakeEnvironment.yaml` (see ***Note***)


```
conda env create -f /<your_path_to>/radseq-preprocessing-pipeline/setup/snakemakeEnvironment.yaml

conda activate radPreprocessing

##check snakemake installed correctly

snakemake --version
```

You may need to close and reopen your terminal after the `conda env create ...` command, or use `source ~/.bashrc` depending on your system.



***Note*** If you already have a snakemake (or suitable Python) installation and would like to avoid installing again, ensure that all of the above software are in your `PATH`.  If you do this instead of installing from the provided radPreprocessing environment (`snakemakeEnvironment.yaml`), you will still need at least the base conda installed/activated - as it's required to handle software dependencies of all the tools used within the workflow itself*
<br>

---

<div id="configure-sample"></div>


## Modifying our configuration files

There are only two files that you as the user are required to modify according to your specific datasets and computing resources available.

These are:

- `config.yaml`
  * Descriptions and details are inside the file.

- `FlexbarSheet.tsv`
  * Mapping of barcodes to their respective replicates, for each library.
  * Provide a `#<Run_name>` at the top to identify results folder


## Running the workflow


From within the main workflow directory (`radseq-preprocessing-pipeline`) where the `Snakefile` is located (make sure your snakemake environment is activated), run with:

```
snakemake --cores # --use-conda
```

where `--cores #` is the maximum number of cores you want to utilise.  If you choose a number of cores that is less than the number of threads defined in your config.yaml for the different tools, then the number of threads will be scaled down to the number of cores defined in the command.

For example:

If you have defined in your config.yaml `bowtie2Threads = 16` and you run snakemake with the command:

```
snakemake --cores 8 --use-conda
```

Your pipeline will be executed with at most 8 threads per job, and so your bowtie2 mapping step will only use 8 threads for each pair of reads, as opposed to 16.

Alternatively if you use a number of cores greater than the defined number of threads:

```
snakemake --cores 32 --use-conda
```

Then each bowtie2 mapping step will use the defined 16 threads, and on top of this, it will execute two bowtie2 mapping's in parallel, since each can use 16 threads and you have allowed for 32 threads (cores) in total.

## Pipeline failing

The most likely reason a pipeline fails to run the whole way through is because no reads could be assigned to one of the barcode-individual combinations as defined in the FlexbarSheet.  Navigate to the results folder and find the ID's for which the fastq files are empty (i.e. file size is `0`).  Remove these ID's from the Flexbarsheet, together with the fulle results folder of the failed run, and re-execute.

## Results
Please look inside your `"/<pathto>/<results>/<Run_name>_0_PreProcessing"` for all results. You will also find a new updated `FlexbarSheetWithCounts.tsv` where you have the number of reads assigned to each barcode combination following inline-barcode demultiplexing.  
