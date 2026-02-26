## Quality Control

Change the parameters based on the organism of interest

## Flagger

```bash
 minimap2 -t 24 --secondary=yes -ax map-pb --MD $Genome_prefix.final.masked.fa .$HiFi_file \
  | samtools view -hSb - \
  | samtools sort -o aln.s.bam


docker run -it --rm -v ${WORKING_DIR}:${WORKING_DIR} mobinasri/flagger:v1.1.0 \
  bam2cov --bam ${WORKING_DIR}/aln.s.bam \
          --output ${WORKING_DIR}/coverage_file.cov.gz \        
          --annotationJson ${WORKING_DIR}/annotations_path.json \
          --threads 16 \                                                 
          --baselineAnnotation whole_genome

mkdir -p ${WORKING_DIR}/hmm_flagger_outputs
docker run -it --rm -v${WORKING_DIR}:${WORKING_DIR} mobinasri/flagger:v1.1.0 \
        hmm_flagger \
            --input ${WORKING_DIR}/coverage_file.cov.gz \
            --outputDir ${WORKING_DIR}/hmm_flagger_outputs  \
            --alphaTsv ${WORKING_DIR}/alpha_optimum_trunc_exp_gaussian_w_4000_n_50.tsv \ ##available in the flagger github page
            --labelNames Err,Dup,Hap,Col \
            --threads 16
```


## BUSCO

```bash
busco -i $Genome_prefix.final.masked.fa -o busco_camelea -l eukaryota -c 8 -m geno -f 
busco -i $Genome_prefix.final.masked.fa -o busco_camelea -l metazoa -c 8 -m geno -f 
busco -i $Genome_prefix.final.masked.fa -o busco_camelea -l mollusca -c 8 -m geno -f 

busco -i $Genome_prefix.hic.a.purged.fasta -o busco_camelea -l eukaryota -c 8 -m geno -f 
busco -i $Genome_prefix.hic.a.purged.fasta -o busco_camelea -l metazoa -c 8 -m geno -f 
busco -i $Genome_prefix.hic.a.purged.fasta -o busco_camelea -l mollusca -c 8 -m geno -f 
```

## merqury

```bash
meryl k=21 count output reads1.meryl DNA_short_1.fq
meryl k=21 count output reads2.meryl DNA_short_2.fq
meryl union-sum output genome.meryl read*.meryl
## run merqury
merqury.sh genome.meryl primary.fasta alternate.fasta final
```
