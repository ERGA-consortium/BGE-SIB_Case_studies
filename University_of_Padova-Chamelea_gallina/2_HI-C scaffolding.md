# HI-C Scaffolding

## HI-C reads mapping

```bash
bwa mem -t 24 $Genome_prefix.hic.p.purged.fa $Hic_file1 | samtools view -@ 10 -Sb - > forward.bam
bwa mem -t 24 $Genome_prefix.hic.p.purged.fa $Hic_file2 | samtools view -@ 10 -Sb - > reverse.bam

```


## Arima preprocessing

```bash
samtools view -h forward.bam | perl /home/mapping_pipeline/filter_five_end.pl | samtools view -Sb - > filter_forward.bam
samtools view -h reverse.bam | perl /home/mapping_pipeline/filter_five_end.pl | samtools view -Sb - > filter_reverse.bam
perl /home/mapping_pipeline/two_read_bam_combiner.pl filter_forward.bam filter_reverse.bam samtools 10 | samtools view -bS -t $Genome_prefix.hic.p.purged.fa.fai - | samtools sort -@ 10 -o sra.bam -
java -Xmx128G -Djava.io.tmpdir=temp/ -jar /home/gatk/picard.jar AddOrReplaceReadGroups INPUT=sra.bam OUTPUT=paired.sra.bam ID=CGAL LB=CGAL SM=CGAL1 PL=ILLUMINA PU=none
java -Xmx128G -XX:-UseGCOverheadLimit -Djava.io.tmpdir=temp/ -jar /home/gatk/picard.jar MarkDuplicates INPUT=paired.sra.bam OUTPUT=$Genome_prefix.primary.bam METRICS_FILE=metrics.cgal.txt TMP_DIR=temp/ ASSUME_SORTED=TRUE VALIDATION_STRINGENCY=LENIENT REMOVE_DUPLICATES=TRUE
samtools index cgal.primary.bam
perl /home/mapping_pipeline/get_stats.pl $Genome_prefix.primary.bam > $Genome_prefix.primary.bam.stats
bamToBed -i $Genome_prefix.primary.bam > alignment.bed
sort -k 4 alignment.bed > tmp && mv tmp alignment.bed
```


## SALSA pipeline

```bash
python /home/SALSA/run_pipeline.py -a $Genome_prefix.hic.p.purged.fa -l $Genome_prefix.hic.p.purged.fa.fai -b alignment.bed -e GATC -o scaffolds -m yes -i 20 -p yes
```

## yagcloser

```bash
minimap2 -t 120 --secondary=yes -ax map-pb scaffolds_FINAL.fasta hifi_reads.fastq | samtools view -hSb - | samtools sort - > aln.s.bam
samtools index aln.s.bam
/home/asset-master/src/detgaps scaffolds_FINAL.fasta > final.gaps.bed
python /home/yagcloser/yagcloser.py -g scaffolds_FINAL.fasta \
    -a aln.s.bam \
    -b final.gaps.bed \
    -o yagcloser_output \
    -s reference_data \
    -mins 2 -f 20 -mcc 2 -prt 0.25 -eft 0.2 -pld 0.2
python /home/yagcloser/scripts/update_assembly_edits_and_breaks.py \
    -i scaffolds_FINAL.fasta \
    -o scaffolds_v2.fasta \
    -e yagcloser_output/reference_data.edits.txt
```



## Hi-C Contact map generation

```bash
samtools faidx scaffolds_v2.fasta
cut -f1,2 scaffolds_v2.fasta.fai > scaffolds_v2.chrom.sizes
bwa index scaffolds_v2.fasta && bwa mem -5SP -t 145 scaffolds_v2.fasta forward.fq backward.fq | samtools view -@ 5 -bhs - > out.bam
samtools view -h out.bam | pairtools parse -c scaffolds_v2.chrom.sizes -o parsed.pairsam
pairtools sort --nproc 8 -o sorted.pairsm parsed.pairsam
pairtools dedup --mark-dups -o deduped.pairsam sorted.pairsam
pairtools select \
'(pair_type == "UU") or (pair_type == "UR") or (pair_type == "RU")' \
-o filtered.pairsam deduped.pairsam
pairtools split --output-sam output.sam filtered.pairsam
samtools view -h output.sam | PretextMap -o map.pretext --sortby length --sortorder descend --mapq 10
```
