use Snakefile_alignment_variantcalling for aligment and variant calling, run with snakemake

filtering:
"bcftools filter --exclude 'MQ0F>0.1' ... | bcftools filter --exclude 'MQ<40' | bcftools filter --exclude 'F_MISSING>0.2' > ..."
"bcftools view --max-alleles 2 --types snps  --apply-filters PASS ... > ..."

proceed with the Snakefile_popgen script, run with snakemake

plot results with sandgoby.R script



used https://github.com/millanek/Dsuite for Dinvestigate

plotted Dinvestigate results with Plot_dinvestigate script in jupyter notebook (R)



f_st:
"vcftools --gzvcf ... --weir-fst-pop ... --weir-fst-pop ... --fst-window-size 10000 --fst-window-step 5000 --out ..."
you need 3 population pairs for PBS

proceed with Snakefile_PBS, run with snakemake

plot the results (similar to dinvestigate plot)


used https://github.com/reslp/phylociraptor for phylogenomic analyses



