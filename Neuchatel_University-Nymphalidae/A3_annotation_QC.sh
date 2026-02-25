#!/bin/bash

### Set directories
INDIR=$MAINDIR/Annotations/C3_braker_reads
OUTDIR=$MAINDIR/Annotations/C4_annotation_QC_reads

# Function to parallelize over samples
annotQC() {
    SAMPLE=$1
    
    ### Select only the longest isoform for each gene
    ### Use script found there: https://bioinformatics.stackexchange.com/questions/595/how-can-longest-isoforms-per-gene-be-extracted-from-a-fasta-file
    longest_iso.py $INDIR/${SAMPLE}/braker.codingseq > $OUTDIR/${SAMPLE}_longest_iso.codingseq

    ### Run busco on the braker protein output (nucleotides)
    source $SOFTWAREDIR/miniconda3/etc/profile.d/conda.sh
    conda activate busco_env

    GENES=$OUTDIR/${SAMPLE}_longest_iso.codingseq
	busco -i $GENES --download_path $MAINDIR/busco_downloads -l lepidoptera_odb10 --offline -o ${SAMPLE}_busco -m transcriptome -c 16

    conda deactivate

    # Count the number of each feature in the GTF
    GTF=$INDIR/${SAMPLE}/braker.gtf
    cat $GTF | grep -v '^#' | cut -f3 | sort | uniq -c > $OUTDIR/${SAMPLE}_number_of_features.txt

    # Get the percentage of mappable RNAseq reads
    samtools flagstats $INDIR/${SAMPLE}/${SAMPLE}_isoseq.bam > $OUTDIR/${SAMPLE}_mapping_stats.txt

    # Get mean gene and transcript length
    awk '$3 == "gene"' $GTF > $INDIR/${SAMPLE}/braker_genes.gtf
    awk '$3 == "transcript"' $GTF > $INDIR/${SAMPLE}/braker_transcripts.gtf
    awk '{ diff = ($4 > $5 ? $4 - $5 : $5 - $4); print $0 "\t" diff }' $INDIR/${SAMPLE}/braker_genes.gtf > $INDIR/${SAMPLE}/braker_genes_length.txt
    awk '{ diff = ($4 > $5 ? $4 - $5 : $5 - $4); print $0 "\t" diff }' $INDIR/${SAMPLE}/braker_transcripts.gtf > $INDIR/${SAMPLE}/braker_transcripts_length.txt
    awk '{ sum += $NF; count++ } END { if (count > 0) print sum / count }' $INDIR/${SAMPLE}/braker_transcripts_length.txt

}

### Parallelize over individuals
export INDIR OUTDIR SAMPLE
export -f annotQC
parallel --colsep '\t' 'annotQC {1}' :::: $MAINDIR/Annotations/C0_rna_QC/samples_rna.txt

### SOFTWARE VERSIONS
# BUSCO 5.7.1
# samtools 1.3.1
