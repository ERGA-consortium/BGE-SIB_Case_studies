#!/bin/bash

# Submit with: bsub -q normal -n 32 -M240000 -J C5_entap -o $MAINDIR/Annotations/C5_entap/C5_entap.out -e $MAINDIR/Annotations/C5_entap/C5_entap.err -R "select[mem>240000] rusage[mem=240000] span[hosts=1]" bash C5_entap.sh

### Set directories
INDIR=$MAINDIR/Annotations/C3_braker_reads
OUTDIR=$MAINDIR/Annotations/C5_entap
DMND1=$SOFTWAREDIR/EnTAP-1.0.1/database/config/bin/eggnog_proteins.dmnd
DMND2=$SOFTWAREDIR/EnTAP-1.0.1/database/config/bin/uniprot_sprot.dmnd
INI=$SOFTWAREDIR/EnTAP-1.0.1/entap_config.ini


# Function to parallelize over samples
entap() {
    SAMPLE=$1
    PRED=$INDIR/${SAMPLE}/braker.aa

    source $SOFTWAREDIR/miniconda3/etc/profile.d/conda.sh
    conda activate braker_env

    mkdir $OUTDIR/${SAMPLE}

    EnTAP --runP -i $PRED -d $DMND1 -d $DMND2 --ini $INI -t 16 --out-dir $OUTDIR/${SAMPLE}
    
    conda deactivate

}

### Parallelize over individuals
export INDIR OUTDIR SAMPLE DMND1 DMND2 INI
export -f entap
parallel --colsep '\t' 'entap {1}' :::: $MAINDIR/Annotations/C0_rna_QC/samples_rna.txt

### SOFTWARE VERSIONS
# EnTAP version 1.0.1
