**# Creates an environment for the WGS pipeline #**





\### **1. Install micromamba** ###



*bash*

*"${SHELL}" <(curl -L micro.mamba.pm/install.sh)*



---------------------------------------------------------------------------------



\### **2. Micromamba for your shell (only once)** ###



&nbsp;	#### **2.1. Initialize** ####



*bash*

*micromamba shell init -s bash -p ~/.micromamba*





	*####* <b>2.2. Restart terminal or run ####</b>



*bash*

*source ~/.bashrc*



---------------------------------------------------------------------------------



\### **3. Create one environment for the pipeline** ###



*bash*

*micromamba create -n genomics \\*

 *-c conda-forge -c bioconda \\*

  *bwa=0.7.19 \\*

  *samtools=1.21 \\*

  *gatk4=4.6.2.0 \\*

  *plink2=2.00a3 \\*

  *vcftools=0.1.17 \\*

  *bcftools=1.21 \\*

  *busco=5.8.2 \\*

  *fastp=0.22.0 \\*

  *fastqc=0.12.1 \\*

  *multiqc=1.17*



---------------------------------------------------------------------------------



\### **4. Activate environment** ###



*bash*

*micromamba activate genomics*



*micromamba install -c bioconda bwa=0.7.19 samtools=1.21 gatk4=4.6.2.0 plink2=2.00a3 vcftools=0.1.17 bcftools=1.21 busco=5.8.2 fastp=0.22.0 fastqc=0.12.1 multiqc=1.17*



---------------------------------------------------------------------------------

