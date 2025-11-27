# Notes: Mind the gaps - Populus tremula pipeline

## Build the environment
Miniforge3-Linux-x86_64 was used as conda/mamba manager

```
# make the environment
mamba create -n chromcomp python=3.11 -y

# activate 
mamba activate chromcomp

# install packages
mamba install -c bioconda \
	bwa=0.7.19 \
	samtools=1.21 \
	gatk4=4.6.2.0 \
	plink2=2.00a3 \
	vcftools=0.1.17 \
	bcftools=1.21 \
	busco=5.8.2 \
	repeatmasker=4.1.8 \
	fastqc=0.12.1 \
	multiqc=1.29 \
	fastp=1.0.1 \
	-y
```


## Soft masking repeats

The masking of repeats was run in two steps
1. First masking using the standard `viridiplantae` database

```
RepeatMasker \
	-pa 28 \
	-species "viridiplantae" \
	-xsmall \
	-gff \
	-dir assemblies/${asm}/ \
	assemblies/${asm}/${asm}.fa
```
After this step add `round-1` to all output files.

2. Using the masked output from the 1st step masking with a species curated repeat database (`Repeats_Aspen_1.0.fna`).  

```
RepeatMasker \
	-pa 28 \
	-engine ncbi \
	-lib ${curated_lib} \
	-xsmall \
	-gff \
	-dir assemblies/${asm}/ \
	assemblies/${asm}/${asm}_round-1_masked.fa
```

## Read quality, filtering and trimming
Read quality assessment
```
fastqc \
    -t 8 \
    reads/raw/*.fq.gz \
    -o qc/fastqc

multiqc \
    qc/fastqc \
    -o qc
```

Quality and adapter trimming
```
fastp \
    -i reads/raw/${sample_id}_R1.fq.gz \
    -I reads/raw/${sample_id}_R2.fq.gz \
    -o reads/trimmed/${sample_id}_R1.trim.fq.gz \
    -O reads/trimmed/${sample_id}_R2.trim.fq.gz \
    -q 20 \
    -u 30 \
    -l 50 \
    -g \
    -h qc/fastp/html/${sample_id}.html \
    -j qc/fastp/json/${sample_id}.json
```

## Read mapping

Indexing
```
# bwa index
bwa index assemblies/${asm}/${asm}_round-2_masked.fa

# .fai index
samtools faidx assemblies/${asm}/${asm}_round-2_masked.fa

# .dict, for GATK.
gatk CreateSequenceDictionary -R assemblies/${asm}/${asm}_round-2_masked.fa
```

Mapping
```
# map bwa mem
bwa mem \
    -M \
    -t 8 \
    -R "@RG\tID:${sample_id}\tSM:${sample_id}\tPL:NovaSeq\tPI:269" \
    assemblies/${asm}/${asm}_round-2_masked.fa \
    reads/trimmed/${sample_id}_R1.trim.fq.gz \
    reads/trimmed/${sample_id}_R2.trim.fq.gz |  \
samtools sort \
    -@8 \
    -o bam/${asm}/${sample_id}.${asm}.bam

# mark duplicates picard
gatk MarkDuplicates \
  -I bam/${asm}/${sample_id}.${asm}.bam \
  -O bam/${asm}/${sample_id}.${asm}.dedup.bam \
  -M qc/dedup/${sample_id}.${asm}.dup.txt

# index
samtools index \
    bam/${asm}/${sample_id}.${asm}.dedup.bam
```

Extract flagstats for the mapping
```
samtools flagstat \
    -O tsv \
    bam/${asm}/${sample_id}.${asm}.bam > qc/flagstat/${sample_id}.${asm}.flagstat.txt
```

## Haplotype calling
Haplotype calling was done on subsets of regions specified by the -L parameter which held a bed file with regions to call.

For the draft assembly, the original 204318 contigs were divided into 20 bed files of ~10220 contigs in each.

For the chromosome assembly, each of the 19 chromosomes were called individually while all scaffolds were combined into one call.

Call per sample and region
```
gatk --java-options "-Xmx7g" HaplotypeCaller \
  -I bam/${asm}/${sample_id}.${asm}.dedup.bam \
  -R assemblies/${asm}/${asm}_round-2_masked.fa \
  --emit-ref-confidence GVCF \
  -L ${bed_file} \
  -O gvcf/${asm}/sample_calls/${sample_id}_${asm}_${run_id}.g.vcf.gz

```

For us the recommended `CombineGVCFs` was very slow to run. Instead we used `GenomicsDBImport` to combine the samples calls per called region.

Combine sample calls

```
gatk --java-options "-Xmx20g -DGATK_STACKTRACE_ON_USER_EXCEPTION=true" GenomicsDBImport \
    --genomicsdb-workspace-path gvcf/${asm}/genomicDB/range_${run_id} \
    --sample-name-map gvcf/${asm}/genomicDB/${asm}_${run_id}_cohort-files.txt \
	--merge-contigs-into-num-partitions 25 \
    --merge-input-intervals \
	--bypass-feature-reader \
	-L ${bed_file}
```

Cohort genotype calling 
```
cd gvcf/${asm}/genomicDB

gatk --java-options "-Xmx20g" GenotypeGVCFs \
    -R ${ref} \
    -V gendb://range_${run_id} \
    -O cohort.${asm}.${run_id}.raw.vcf.gz
```

Concat all regions called files into one vcf file

```
bcftools concat \
    --threads 4 \
    --file-list gvcf/${asm}/genomicDB/list_vcf_to_concat.txt \
    -Oz -o gvcf/${asm}/cohort.${asm}.raw.vcf.gz
```

Adding a step to write the index with tabix to get `.tbi` index because bcftools writes a `.csi` index which doesn't work with `gatk VariantFiltration` 

```
tabix \
    -p vcf \
    gvcf/${asm}/cohort.${asm}.raw.vcf.gz
```


Filtering calls
Add a hard filter

```
gatk VariantFiltration \
    -R assemblies/${asm}/${asm}_round-2_masked.fa \
    -V gvcf/${asm}/cohort.${asm}.raw.vcf.gz \
    --filter-expression "QD < 2.0 || FS > 60.0 || MQ < 40.0" \
    --filter-name "basic_snp_filter" \
    -O gvcf/${asm}/cohort.${asm}.flt.vcf.gz
```

Final filtering and retaining only:
- the encoded filter PASS (previous step)
- biallelic SNPs 
- MAF > 0.05 
- maximum of 80% missing data

```
bcftools view \
    -f PASS \
    -m2 \
    -M2 \
    -v snps \
    gvcf/${asm}/cohort.${asm}.flt.vcf.gz | \
vcftools --vcf - \
    --max-missing 0.8 \
    --maf 0.05 \
    --min-alleles 2 \
    --max-alleles 2 \
    --recode \
    --recode-INFO-all \
    --stdout | \
bcftools view -Oz \
    -o gvcf/${asm}/cohort.${asm}.final.recode.vcf.gz

tabix \
    -p vcf \
    gvcf/${asm}/cohort.${asm}.final.recode.vcf.gz
```

## Popgen stats
### PCA
Plink complained about too many contigs so the vcf files were modified to avoid this complaint.
What was done:
1. For each bed file used during the haplotype calling the chromosome was replaced by the bed file number (1-20)
2. To avoid duplicate positions, positions were recalculated as $1$2, where:
    - $1 = original contig number (remove all letters and 0 at start. eg. Potra000019 -> 19) 
    - $2 the original position
    - e.g if the chromosome was Potra000019 it became 19 and with position 1054 the final position would be 191054
3. IDs for each SNP was premade as original_contig/chromosome:new_contig:new_position  eg. contig: Potra000019, position: 1054 -> 19:19:1054. The original contig/chromosome is added because only using new_contig:new_position would create duplicates due to the many contigs
LD prune the data
```
plink2 \
    --vcf gvcf/${asm}/${asm}_plinkMOD.vcf.gz \
    --allow-extra-chr \
    --bad-ld \
   	--indep-pairwise 50 10 0.2 \
  	--out stats/${asm}/${asm}
```

PCA needs a .freq file, make one
```
plink2 \
    --vcf gvcf/${asm}/${asm}_plinkMOD.vcf.gz \
    --allow-extra-chr \
    --extract stats/${asm}/${asm}.prune.in \
    --freq \
  	--out stats/${asm}/freq_${asm}
```

Run PCA
```
plink2 \
    --vcf gvcf/${asm}/${asm}_plinkMOD.vcf.gz \
    --allow-extra-chr \
    --extract stats/${asm}/${asm}.prune.in \
    --read-freq stats/${asm}/freq_${asm}.afreq \
   	--pca 10 \
   	--out stats/${asm}/pca_${asm}
```

### vcftools stats
Nucleotide diversity
```
vcftools \
    --gzvcf gvcf/${asm}/${asm}_final.vcf.gz \
    --site-pi \
    --out stats/${asm}/pi_${asm}
```

Observed heterozygosity
```
vcftools \
    --gzvcf gvcf/${asm}/${asm}_final.vcf.gz \
    --het \
    --out stats/${asm}/het_${asm}
```

Pairwise FST
```
vcftools \
    --gzvcf gvcf/${asm}/${asm}_final.vcf.gz \
    --weir-fst-pop pop1.txt \
    --weir-fst-pop pop2.txt \
    --out stats/${asm}/fst_${asm}
```
pop1.txt (East, higher elevation)
- NordAsp_E17
- NordAsp_E21
- NordAsp_E28
- NordAsp_X701
- NordAsp_X706

pop2.txt (West, lower elevation)
- NordAsp_G17
- NordAsp_M14
- NordAsp_M21
- NordAsp_M29
- NordAsp_M30

### Addition
Pairwise Fst was also calculated in windows for three populations using the same code as above with the added `--fst-window-size 10000`

pop1.txt (West)
- NordAsp_M14
- NordAsp_M21
- NordAsp_M29
- NordAsp_M30
- NordAsp_M48

pop2.txt (Mountain)
- NordAsp_H5
- NordAsp_N8
- NordAsp_N29
- NordAsp_B7

pop3.txt (East)
- NordAsp_L54
- NordAsp_L49
- NordAsp_L46
- NordAsp_F44
- NordAsp_F31
