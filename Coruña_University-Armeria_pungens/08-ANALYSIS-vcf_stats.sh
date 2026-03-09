# CONFIGURATION:

INPUT_FILE="/route/to/vcf/file"
OUTPUT_FILE="/route/to/output/file"  # base name for all outputs

# EXECUTION:

vcftools --vcf "$INPUT_FILE" --site-mean-depth --out "$OUTPUT_FILE"  # depth coverage
vcftools --vcf "$INPUT_FILE" --freq --out "$OUTPUT_FILE"  # alleles per site
vcftools --vcf "$INPUT_FILE" --keep-only-indels --out "$OUTPUT_FILE"  # indels
vcftools --vcf "$INPUT_FILE" --extract-FORMAT-info GQ --out "$OUTPUT_FILE"  # quality (GQ)
vcftools --vcf "$INPUT_FILE" --missing-site --out "$OUTPUT_FILE"  # missing data (SNP)
vcftools --vcf "$INPUT_FILE" --missing-indv --out "$OUTPUT_FILE"  # missing data (sample)