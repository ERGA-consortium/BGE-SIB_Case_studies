#!/bin/bash -l

#SBATCH -A hpc2n2025-120	
#SBATCH -n 1
#SBATCH -c 1
#SBATCH -t 1-00:00:00
#SBATCH -J pca
#SBATCH --output=/proj/nobackup/hpc2nstor2025-059/mimmi/mind_the_gaps/reports/sbatch/pca/sbatch_R-%x_%j-%a.out
#SBATCH --error=/proj/nobackup/hpc2nstor2025-059/mimmi/mind_the_gaps/reports/sbatch/pca/sbatch_R-%x_%j-%a.err
#SBATCH -a 0-1

eval "$(mamba shell hook --shell bash)"
mamba activate chromcomp

cd /proj/nobackup/hpc2nstor2025-059/mimmi/mind_the_gaps

run_list=(draft chrom)
asm=$(echo ${run_list[$SLURM_ARRAY_TASK_ID]})

species=$1

echo ${species}
echo ${asm}

# modify chromosome numbers for plink so it does not complain of too many contigs
# for each bedfile (20 in total) replace chromosome with bedfile number. 
# for position (column 2), remove all letters and 0 at the start of the chrom then make the position a combination of chrom number and position ($1$2)
# for ID (column 3) combine old_chromosome/contig":"new_contig":"newposition ($1$2)
echo "modifying vcf for plink"

mkdir data/${species}/gvcf/${asm}/temp

bcftools view \
    --no-version -h \
    data/${species}/gvcf/${asm}/cohort.${asm}.final.recoded.vcf.gz | \
    sed '/^##contig/d'> data/${species}/gvcf/${asm}/temp/header.vcf
tail -n1 data/${species}/gvcf/${asm}/temp/header.vcf > data/${species}/gvcf/${asm}/temp/colnames.vcf
sed '/^#CHROM/d' -i data/${species}/gvcf/${asm}/temp/header.vcf

for bed_file in $(ls data/${species}/bedfiles/${asm}/*.bed | sort -V)
do 
    echo ${bed_file}
    bed_id=$(echo $(basename ${bed_file} | sed 's/[a-zA-Z_.]//g'))
    bcftools view \
        --no-version -H \
        -R ${bed_file} \
        data/${species}/gvcf/${asm}/cohort.${asm}.final.recoded.vcf.gz | \
    cut -f1,2 | \
    tail -n1 | \
    sed 's/^[a-zA-Z]*//' | \
    sed "s/^0*/${bed_id}\t/" | \
    awk '{print "##contig=<ID="$1",length="$1$3+10">"}' >> data/${species}/gvcf/${asm}/temp/header.vcf

    bcftools view \
        --no-version -H \
        -R ${bed_file} \
        data/${species}/gvcf/${asm}/cohort.${asm}.final.recoded.vcf.gz | \
    sed 's/^[a-zA-Z]*//' | \
    sed 's/^0*//' | \
    awk -v b="${bed_id}" '{print b,b$2,$1":"b":"b$2,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20}' OFS='\t' >> data/${species}/gvcf/${asm}/temp/variants.vcf
done

cat data/${species}/gvcf/${asm}/temp/header.vcf \
    data/${species}/gvcf/${asm}/temp/colnames.vcf \
    data/${species}/gvcf/${asm}/temp/variants.vcf | \
bcftools sort | \
bcftools view \
    --no-version \
    -Oz -o data/${species}/gvcf/${asm}/${asm}_plinkMOD.vcf.gz \

tabix --csi -p vcf data/${species}/gvcf/${asm}/${asm}_plinkMOD.vcf.gz

rm -rf data/${species}/gvcf/${asm}/temp

# Principal component analysis
# Step A – LD‑prune at r² < 0.2, 50 kb windows, sliding 10 variants
echo "LD pruning"
plink2 --vcf data/${species}/gvcf/${asm}/${asm}_plinkMOD.vcf.gz \
    --allow-extra-chr \
    --bad-ld \
   	--indep-pairwise 50 10 0.2 \
  	--out data/${species}/stats/${asm}/${asm}

# Step B make freq file
echo "create a .freq file"
plink2 --vcf data/${species}/gvcf/${asm}/${asm}_plinkMOD.vcf.gz \
    --allow-extra-chr \
    --extract data/${species}/stats/${asm}/${asm}.prune.in \
    --freq \
  	--out data/${species}/stats/${asm}/freq_${asm}

# Principal component analysis 
# Step C - PCA on the pruned SNP subset
echo "run plink pca"
plink2 --vcf data/${species}/gvcf/${asm}/${asm}_plinkMOD.vcf.gz \
    --allow-extra-chr \
    --extract data/${species}/stats/${asm}/${asm}.prune.in \
    --read-freq data/${species}/stats/${asm}/freq_${asm}.afreq \
   	--pca 10 \
   	--out data/${species}/stats/${asm}/pca_${asm}

rm data/${species}/gvcf/${asm}/${asm}_plinkMOD.vcf.gz*