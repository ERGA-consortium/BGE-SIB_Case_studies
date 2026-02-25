#!/bin/bash

### Set directories
INDIR=$MAINDIR/assemblies/chromosomes_only
OUTDIR=$MAINDIR/Annotations/C2_repeat_annotation

### RepeatMasker in parallel
rmasker() {
    SAMPLE=$1
    REF=$INDIR/${SAMPLE}.fasta
    # Erebia repeat library previously generated using EarlGrey
    LIB=$MAINDIR/Annotations/C1b_blast_lib/Erebia_lib_clstrd_hostgenesfilt.fasta

    source $SOFTWAREDIR/miniconda3/etc/profile.d/conda.sh
    conda activate earlgrey_env

    mkdir $OUTDIR/${SAMPLE}
    cd $OUTDIR/${SAMPLE}
    RepeatMasker -pa 6 -e rmblast -lib $LIB -dir $OUTDIR/${SAMPLE} -xsmall -html -gff $REF
    # -pa runs in parallel (pa 4 uses 16 cores in total in combination with rmblast)
    # -e specifies the search engine (rmblast = ncbi is the one used by EarlGrey)
    # -xsmall for softmasking (repeats in lowercase) rather than changing them to Ns
    # -html and -gff create additional output files

    conda deactivate

}

### Parallelize over individuals
export INDIR OUTDIR SAMPLE
export -f rmasker
parallel --colsep '\t' 'rmasker {1}' :::: $MAINDIR/After_assemblies/samples_codes_full.txt

### SOFTWARE VERSIONS
# RepeatMasker version 4.1.5
