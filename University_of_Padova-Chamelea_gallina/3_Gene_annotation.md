# Gene annotation

## hisat2

```bash
hisat2-build -p 10 $Genome_prefix.final.masked.fa $Genome_prefix 
hisat2 -p 20 -x $Genome_prefix -1 Gill_1.fq.gz -2 Gill_2.fq.gz | samtools sort -@ 10 > Gill.bam &
hisat2 -p 20 -x $Genome_prefix -1 Gonad_1.fq.gz -2 Gonad_2.fq.gz | samtools sort -@ 10 > Gonad.bam &
hisat2 -p 20 -x $Genome_prefix -1 Mantle_1.fq.gz -2 Mantle_2.fq.gz | samtools sort -@ 10 > Mantle.bam &
hisat2 -p 20 -x $Genome_prefix -1 Haemolymph_1.fq.gz -2 Haemolymph_2.fq.gz | samtools sort -@ 10 > Haemolymph.bam 
hisat2 -p 20 -x $Genome_prefix -1 Hepatopancreas_1.fq.gz -2 Hepatopancreas_2.fq.gz | samtools sort -@ 10 > Hepatopancreas.bam 
```

## Braker3

```bash
singularity exec singularity/braker3.sif braker.pl --genome $Genome_prefix.final.masked.fa \
--prot_seq protein.fa \
--bam Gill.bam,Gonad.bam,Mantle.bam,Haemolymph.bam,Hepatopancreas.bam \
--threads 18  --softmasking --species={$Genome_prefix} --skip_fixing_broken_genes 
```

## Trinity 
```bash
samtools merge -o all.bam Gill.bam Gonad.bam Mantle.bam Haemolymph.bam Hepatopancreas.bam 
samtools sort -o all.sorted.bam all.bam -@ 12
Trinity --genome_guided_bam all.sorted.bam --max_memory 250G --genome_guided_max_intron 20000 --CPU 96
```

## PASA
```bash
Launch_PASA_pipeline.pl -c pasa.alignAssembly.Template.txt -R -g $Genome_prefix.final.masked.fa -t Trinity-GG.fasta -C -r --ALIGNERS gmap,minimap2,blat
```

## Evidence Modeler
```bash
touch weights.txt
echo -e "ABINITIO_PREDICTION\tAUGUSTUS\t1" >> weights.txt
echo -e "TRANSCRIPT\tassembler-mydb.sqlite\t10" >> weights.txt
EVidenceModeler --CPU 32 --sample_id cgal --weights weights.txt --genome $Genome_prefix.final.masked.fa --gene_predictions augustus.evm.gff3 --transcript_alignments pasa_db.pasa_assemblies.gff3 --segmentSize 100000 --overlapSize 10000
```

## Add UTRs with PASA (2 iterations)
```bash
Load_Current_Gene_Annotations.dbi -c pasa.alignAssembly.Template.txt -g /home/umberto/raid5/camelea/softmasked_xs/Cgal.final.fa.masked -P cgal.EVM.gff3
PASApipeline/misc_utilities/pasa_gff3_validator.pl cgal.EVM.gff3
Launch_PASA_pipeline.pl -c annotCompare.config -A -g $Genome_prefix.final.masked.fa  -t Trinity-GG.fasta
```
## InterProScan
```bash
./interproscan.sh -cpu 36 -d interproscan -dp -goterms -iprlookup -i cgal.EVM.pep
```
## blastp

```bash
blastp \
-query cgal.EVM.pep \
-db swissprot \
-out cgal.EVM.pep.swissprot.txt \
-outfmt "10 evalue length qseqid qlen qstart qend sacc slen sstart send pident nident sstrand qcovs qseq sseq sgi stitle" \
-num_threads 32 -evalue 1e-5 -max_target_seqs 1 \
-max_hsps 1
```
## Non coding RNA
```bash
cmpress Rfam.cm
esl-seqstat $Genome_prefix.final.masked.fa
cmscan -Z 3615.7 --cut_ga --rfam --nohmmonly --tblout mrum-genome.tblout --fmt 2 --clanin Rfam.clanin Rfam.cm $Genome_prefix.final.masked.fa >$Genome_prefix.genome.cmscan
```
