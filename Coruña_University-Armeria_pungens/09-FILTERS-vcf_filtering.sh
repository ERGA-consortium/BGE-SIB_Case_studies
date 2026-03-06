# The script is designed to be executed step by step!!!

# CONFIGURATION:

INPUT_DIR="/route/to/raw/vcf/file/obtained/with/stacks/directory"
FILTER_DIR="/route/to/filtering/directory"
BCFTOOLS_DIR="/route/to/bcftools-1.22/program/files/directory"
NO_PARALOG_FILE="/route/to/no/paralog/snps/file"
MONO_FILE="/route/to/monomorphic/snps/file"  # see below how to calculate it with R 

# EXECUTION:

# 1) Minimum depth coverage filter. Values: 7 or 10

# CAUTION! Sorting is executed with bcftools in PATH, but FORMAT/DP is executed with bcftools executable file

bgzip -c "${INPUT_DIR}/populations.snps.vcf" > "${FILTER_DIR}/snps.vcf.gz"
bcftools sort "${FILTER_DIR}/snps.vcf.gz" -o "${FILTER_DIR}/snps.vcf.gz"  # overwritting
tabix --csi -p vcf "${FILTER_DIR}/snps.vcf.gz"
export BCFTOOLS_PLUGINS="${BCFTOOLS_DIR}/plugins"
"${BCFTOOLS_DIR}/bcftools" +setGT "${FILTER_DIR}/snps.vcf.gz" -O v -o "${FILTER_DIR}/snps_dp.vcf" -- -t q -n. -i 'FORMAT/DP < 7'

# 2) Minimum quality filter. Value: 30

vcftools --vcf "${FILTER_DIR}/snps_dp.vcf" --minGQ 30 --recode --recode-INFO-all --out "${FILTER_DIR}/snps_dp_gq"

# 3) Paralog filter. Value: 0.6 (H max) and 15 (|D|), depending on HDplot (see https://github.com/gjmckinney/HDplot)

awk -F ',' 'NR > 1 && $2+0 > 0 {print $1 "\t" $2}' "$NO_PARALOG_FILE" > "${FILTER_DIR}/no_paralog_snps_list"
vcftools --vcf "${FILTER_DIR}/snps_dp_gq.recode.vcf" --positions "${FILTER_DIR}/no_paralog_snps_list" --recode --recode-INFO-all --out "${FILTER_DIR}/snps_dp_gq_noparalog"

# 4) Hardy-Weinberg equilibrium deviation filter (this filter can be omitted: rename inputs/outputs in consequence!). Values: - (omit), 0.01 or 0.05

vcftools --vcf "${FILTER_DIR}/snps_dp_gq_noparalog.recode.vcf" --hwe 0.05 --recode --recode-INFO-all --out "${FILTER_DIR}/snps_dp_gq_noparalog_hwe"

# 5) Minimum allele frequency filter (this filter can be omitted: rename inputs/outputs in consequence!). Values: - (omit this filter), 0.01, 0.02 or 0.05

vcftools --vcf "${FILTER_DIR}/snps_dp_gq_noparalog_hwe.recode.vcf" --maf 0.01 --recode --recode-INFO-all --out "${FILTER_DIR}/snps_dp_gq_noparalog_hwe_maf"

# 6) Maximum missing data (SNP) filter. Values: 0.95, 0.85, 0.75 or 0.5

vcftools --vcf "${FILTER_DIR}/snps_dp_gq_noparalog_hwe_maf.recode.vcf" --max-missing 0.5 --recode --recode-INFO-all --out "${FILTER_DIR}/snps_dp_gq_noparalog_hwe_maf_misssnp"

# 7) Maximum missing data (individual) filter. Values: 0.05, 0.15, 0.25 or 0.5

vcftools --vcf "${FILTER_DIR}/snps_dp_gq_noparalog_hwe_maf_misssnp.recode.vcf" --missing-indv --out "${FILTER_DIR}/miss_ind"
cat "${FILTER_DIR}/miss_ind.imiss" | awk '$5+0 > 0.5' | cut -f 1 > "${FILTER_DIR}/miss_ind_list.txt"
vcftools --vcf  "${FILTER_DIR}/snps_dp_gq_noparalog_hwe_maf_misssnp.recode.vcf" --remove "${FILTER_DIR}/miss_ind_list.txt" --recode --recode-INFO-all --out "${FILTER_DIR}/snps_dp_gq_noparalog_hwe_maf_misssnp_missind"

# 8) Monomorphic loci filter

# Identifying monomorphic SNPs (with R):

# library(gaston,quietly=TRUE)
# vcfInput<-read.VCF("/route/to/vcf/file")
# mono_snps <- vcfInput@snps[vcfInput@snps$maf == 0, ]
# name = paste0(gsub("\\.vcf$", "", "/route/to/vcf/file"), "_mono.txt")
# write.table(mono_snps$id, file = name, quote = F, row.names = F, col.names = F)

vcftools --vcf  "${FILTER_DIR}/snps_dp_gq_noparalog_hwe_maf_misssnp_missind.recode.vcf" --exclude "$MONO_FILE" --recode --recode-INFO-all --out "${FILTER_DIR}/snps_dp_gq_noparalog_hwe_maf_misssnp_missind_nomono"
