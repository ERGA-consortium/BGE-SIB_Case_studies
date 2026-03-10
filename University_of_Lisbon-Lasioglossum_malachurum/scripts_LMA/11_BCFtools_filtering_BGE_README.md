\# **SNP filtering using BCFtools on a cohort VCF** #



\\#### 1. Before Running ####



You must edit the script to set your own input/output directories



\\\*vcffile\\\* = \\\*<PATH/to/input\\\_vcf>\\\*



\\\*ref\\\_genome\\\* = \\\*<PATH/to/reference\\\_genome>\\\*



\\\*outtxt\\\* = \\\*<output\\\_text\\\_file\\\_name.txt\\>\*



---------------------------------------------------------------------------------



\\#### 2. Parameters used for bwa-mem ####



\\\*-i TYPE="snp"\\\* (keep only SNPs)



\\\*--fasta-ref <ref>\\\* (reference genome)



\\\*-m -snps\\\* (multiallelic calling mode)



\\\*--min-alleles 2 --max-alleles 2\\\* (keeps snps with exaclty two alleles)



\\\*-i 'QUAL>=20'\\\* (keep snps with QC>=20)



\\\*-e 'F\_MISSING>0.8' 'F\_MISSING>0.5'\\\* (exlude SNPs based on missing data)



\\\*-e 'GT\[\*] = "miss"'\\\* (remove SNPs containing genotype states equal to "miss")



---------------------------------------------------------------------------------



\\#### \\\*\\\*3. Activate the micromamba environment containing GATK4 BCFtools\\\*\\\* ####



\\\*bash\\\*

\\\*micromamba activate <your\\\_environment\\\_name>\\\*



---------------------------------------------------------------------------------



\\#### 4. Run the script ####



\\\*bash 11\\\_BCFtools\\\_filtering\\\_BGE.sh\\\*



---



