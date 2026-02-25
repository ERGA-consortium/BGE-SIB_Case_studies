#!/bin/bash

### Set directories
INDIR=$MAINDIR/Annotations/C2_repeat_annotation
OUTDIR=$MAINDIR/Annotations/C3b_braker_reads

# Function to parallelize over samples
braker() {
    SAMPLE=$1
    CODE=$2
    SAMPLEREF=$3
    PROT=$OUTDIR/Arthropoda.fa #odb11
    READS=$MAINDIR/raw_data/RNAseq/concat/${SAMPLE}_reads.fastq
    
    ### Copying reference genome in wd otherwise Braker does not see it
    mkdir $OUTDIR/${SAMPLE}
    cp $INDIR/${SAMPLEREF}/${SAMPLEREF}.fasta.masked $OUTDIR/${SAMPLE}
    REF=$OUTDIR/${SAMPLE}/${SAMPLEREF}.fasta.masked

    ### Map the reads to the reference
    minimap2 -t 16 -ax splice:hq -uf $REF $READS > $OUTDIR/${SAMPLE}/isoseq.sam     
    samtools view -bS --threads 16 $OUTDIR/${SAMPLE}/isoseq.sam -o $OUTDIR/${SAMPLE}/${SAMPLE}_isoseq.bam
    rm $OUTDIR/${SAMPLE}/isoseq.sam

    ### Pull the container (only once)
    # singularity build braker3_lr.sif docker://teambraker/braker3:isoseq
    
    ### Run BRAKER adapter for long reads + protein db
    singularity exec -B ${PWD}:${PWD} $SOFTWAREDIR/braker3_lr.sif braker.pl --workingdir=$OUTDIR/${SAMPLE} \
        --genome=$REF --prot_seq=$PROT --useexisting --species=${SAMPLE}_2 --bam=$OUTDIR/${SAMPLE}/${SAMPLE}_isoseq.bam --threads=16

}

### Parallelize over individuals
export INDIR OUTDIR SAMPLE CODE SAMPLEREF
export -f braker
parallel --colsep '\t' 'braker {1} {2} {3}' :::: $MAINDIR/Annotations/C0_rna_QC/samples_rna.txt

### SOFTWARE VERSIONS
# minimap2 v2.21
# BRAKER version 3.0.8 (for long reads + prot db)
