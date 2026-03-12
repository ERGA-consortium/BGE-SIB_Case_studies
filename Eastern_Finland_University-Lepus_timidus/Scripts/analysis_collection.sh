#Sex determination
module load samtools/1.21

for i in ../04_bwa_align/hybrids/*bam; do echo $i; samtools view $i NC_084851.1 | wc -l; done > Y.readcounts.txt
cat Y.readcounts.txt |tr '\n' '\t' | sed 's/..\/04_bwa_align\/hybrids\//\n/g' | sed 's/.bam//' > Y.readcounts.list

for i in ../04_bwa_align/hybrids/*bam; do echo $i; samtools view $i NC_084850.1 | wc -l; done > X.readcounts.txt
cat X.readcounts.txt |tr '\n' '\t' | sed 's/..\/04_bwa_align\/hybrids\//\n/g' | sed 's/.bam//' > X.readcounts.list


#Plink
module load plink/2.00a5

plink2 --allow-extra-chr \
  --allow-no-sex \
  --make-pgen \
  --out radseq \
  --vcf radseq.all.filtered.vcf.gz

#Missingness and allele frequency stats
plink2  --allow-extra-chr \
  --freq \
  --geno-counts \
  --genotyping-rate \
  --hardy \
  --missing \
  --out radseq.stats \
  --pfile radseq


#PCA
plink2 --allow-extra-chr \
       --chr-set 26 \
       --pfile radseq \
       --out ragseq_pca \
       --pca 20

#ADMIXTURE
module load bcftools/1.16
module load admixture/1.3.0

#Rename chromosomes to numeric
bcftools annotate --rename-chrs chr.main.txt -o all.chr_rename.vcf.gz all.merged.filtered.vcf.gz
bcftools norm -m-any -o all.merged.chr_rename.biall.vcf.gz all.merged.chr_rename.vcf.gz
#Include only autosomes:
zcat all.merged.chr_rename.biall.vcf.gz | grep -v NW > all.primary.vcf
#convert to plink files while removing samples with high missingness + calculating some stats
plink2 --vcf all.primary.vcf --chr-set 26 --make-bfile --out all.bi --remove ../rem.txt --missing --freq
#Run - repeated for K=2 and K=4 
admixture -j4 -B all.bi.bed 2


#ELAI
module load plink/1.90 #makes conversion to bimbam easier, than with plink2
#Input prep - chrX
for j in LE_FIN LT_FIN LT_parental LE_parental LE_all LT_all; do
    echo $j;
    bcftools view -S ../POPS/nonmix.$j.txt -o X.$j.bi.ann.snp.vcf.gz X.bi.ann.snp.vcf.gz;
done
for j in LE.LT LT.LE LE.any LT.any; do
    echo $j;
    bcftools view -S ../POPS/admix.$j.txt -o admix.X.$j.bi.ann.snp.vcf.gz X.bi.ann.snp.vcf.gz;
done
for j in LT_parental LE.LT LT.LE LE.any LE_all LT_all; do
    echo $j;
    plink --make-bed \
	  --chr-set 26 \
	  --allow-no-sex \
	  --allow-extra-chr \
	  --out X.$j \
	  --double-id \
	  --vcf X.$j.bi.ann.snp.vcf.gz >> perpop.convert_to_plink.log;
done
for j in LE_FIN LT_FIN LE_parental LT_parental LT.LE LE.LT LE_all LT_all; do
    echo $j;
    plink \
	--recode-bimbam \
	--out X.$j \
	--bfile X.$j \
	--snps-only \
	--allow-extra-chr >> bimbam_convert.log;
done

#Split and adjust vcf files of autosomes
for i in {1..23}; do echo $i; bcftools view -r $i -o $i.m05.s075.vcf.gz ALL.m0.5.s0.75.vcf.gz; done
for i in {1..23}; do echo $i; bcftools norm -m -any -o $i.bi.m05.s075.vcf.gz $i.m05.s075.vcf.gz ; done
for i in {1..23}; do echo $i; bcftools annotate -o $i.bi.ann.m05.s075.vcf.gz --set-id +'%CHROM\_%POS\_%REF\_%FIRST_ALT' $i.bi.m05.s075.vcf.gz; done
for i in {1..23}; do echo $i; bcftools view -v snps -o $i.bi.ann.snp.vcf.gz $i.bi.ann.m05.s075.vcf.gz; done

#Non-admixed populations:
for j in LE_FIN LT_FIN LT_parental LE_parental; do for i in {1..23}; do echo $i $j; bcftools view -S ../POPS/nonmix.$j.txt -o $i.$j.bi.ann.snp.vcf.gz $i.bi.ann.snp.vcf.gz; done; done
for j in LT_FIN LE_FIN LE_parental LT_parental; do for i in {1..23}; do echo $i $j; plink --make-bed --chr-set 26 --allow-no-sex --allow-extra-chr --out $i.$j --double-id --vcf $i.$j.bi.ann.snp.vcf.gz >> perpop.convert_to_plink.log; done; done
for j in LE.LT LT.LE LE_all LT_all; do for i in {1..23}; do echo $i $j; plink --make-bed --chr-set 26 --allow-no-sex --allow-extra-chr --out $i.$j --double-id --vcf $i.$j.bi.ann.snp.vcf.gz >> convert_to_plink.$j.log; done; done
for j in LE_FIN LT_FIN LE_parental LT_parental; do for i in {1..23}; do echo $i $j; plink --recode-bimbam --out $i.$j --bfile $i.$j --snps-only; done; done

#Admixed populations:
for i in {1..23}; do echo $i; bcftools view -S ../POPS/admix.LE.any.txt -o admix.$i.LE.bi.ann.snp.vcf.gz $i.bi.ann.snp.vcf.gz; done
for i in {1..23}; do echo $i; plink --make-bed --chr-set 26 --allow-no-sex --allow-extra-chr --out $i --double-id --vcf $i.bi.ann.snp.vcf.gz >> convert_to_plink.log; done
for i in {1..23}; do echo $i; plink --make-bed --chr-set 26 --allow-no-sex --allow-extra-chr --out admix.$i.LE.any --double-id --vcf admix.$i.LE.bi.ann.snp.vcf.gz >> admix.LE.any.perpop.convert_to_plink.log; done
for i in {1..23}; do echo $i; plink --recode-bimbam --out admix.$i.LE.any --bfile admix.$i.LE.any --snps-only; done > LE_any.bimbam.log

#Run example
#ELAI version 1.01.
export PATH="/projappl/project_2002674/programs/:$PATH"
elai-lin -g ../PER_CHROM/1.LT_FIN.recode.geno.txt -p 10 \
	 -g ../PER_CHROM/1.LT_parental.recode.geno.txt -p 11 \
	 -g ../PER_CHROM/1.LE_parental.recode.geno.txt -p 12 -g \
	 ../PER_CHROM/admix.1.LE.any.recode.geno.txt -p 1 \
	 -pos ../PER_CHROM/1.recode.pos.txt \
	 -o 1.LE_any \
	 -C 3 -c 73 -s 30 -sem 0 -mg 1
#Repeated for all chromosomes, and did some reruns with slightly modified parameters



#Vcftools analyses
module load vcftools/0.1.17

#Fst
for i in AMU Hir IN Joe KAM KOL Lap MAG Oulu Pri Pun SCH URA YAK; do
    for j in AMU Hir IN Joe KAM KOL Lap MAG Oulu Pri Pun SCH URA YAK; do
	echo $i $j;
	vcftools --weir-fst-pop $i.pop.list --weir-fst-pop $j.pop.list --fst-window-size 500000 --fst-window-step 1000000 --out $i.$j --vcf LT.filt.vcf ;
    done;
done 2> LT_fst.log > LT_fst.out

for i in ALM Hel Joe Kor LAS Mjo Oulu PYR Pun STY Val; do
    for j in ALM Hel Joe Kor LAS Mjo Oulu PYR Pun STY Val; do
	vcftools --weir-fst-pop $i.pop.list --weir-fst-pop $j.pop.list --fst-window-size 500000 --fst-window-step 1000000 --out $i.$j --vcf LE.filt.vcf ;
    done;
done 2> LE_fst.log > LE_fst.out


for i in AMU Hir IN Joe KAM KOL Lap MAG Oulu Pri Pun SCH URA YAK; do
    for j in AMU Hir IN Joe KAM KOL Lap MAG Oulu Pri Pun SCH URA YAK; do
	echo $i $j;
	vcftools --weir-fst-pop LT_${i}.list --weir-fst-pop LT_${j}.list --fst-window-size 500000 --fst-window-step 100000 --out LT_LT_FST/LT_${i}.LT_${j} --vcf LT.filt.vcf ;
    done;
done 2> LT_LT_FST/LT_fst.log > LT_LT_FST/LT_fst.out

for i in ALM Hel Joe Kor LAS Mjo Oulu PYR Pun STY Val; do
    for j in ALM Hel Joe Kor LAS Mjo Oulu PYR Pun STY Val; do
	vcftools --weir-fst-pop LE_${i}.list --weir-fst-pop LE_${j}.list --fst-window-size 500000 --fst-window-step 100000 --out LE_LE_FST/LE_${i}.LE_${j} --vcf LE.filt.vcf ;
    done;
done 2> LE_LE_FST/LE_fst.log > LE_LE_FST/LE_fst.out

for i in ALM Hel Joe Kor LAS Mjo Oulu PYR Pun STY Val; do
    for j in AMU Hir IN Joe KAM KOL Lap MAG Oulu Pri Pun SCH URA YAK; do
	vcftools --weir-fst-pop LE_${i}.list --weir-fst-pop LT_${j}.list --fst-window-size 500000 --fst-window-step 100000 --out LT_LE_FST/LE_${i}.LT_${j} --gzvcf pass.LTLE_autosomes.miss0.75.maf0.05.samp0.25.LepEur.recode.vcf.gz ;
    done;
done 2> LT_LE_FST/LELT.fst.log > LT_LE_FST/LELT.fst.out

for j in AMU Hir IN Joe KAM KOL Lap MAG Oulu Pri Pun SCH URA YAK; do
    for i in ALM Hel Joe Kor LAS Mjo Oulu PYR Pun STY Val; do
	vcftools --weir-fst-pop LISTS_MAPS/LT_${j}.list --weir-fst-pop LISTS_MAPS/LE_${i}.list --fst-window-size 500000 --fst-window-step 100000 --out LT_LE_FST/LT_${j}.LE_${i} --gzvcf pass.LTLE_autosomes.miss0.75.maf0.05.samp0.25.LepEur.recode.vcf.gz ;
    done;
done 2> LT_LE_FST/LTLE.fst.log > LT_LE_FST/LTLE.fst.out

#Summarizing results - example command
grep -E "\-\-out | estimate" LELT.fst.log | grep -v Outputting | sed 's/\//:/' | cut -d':' -f2- | paste - - - | sed 's/\./\t/' > LELT_subpop_fst_summary.tsv

#Tajima-D - same populations and window sizes, with --TajimaD option


