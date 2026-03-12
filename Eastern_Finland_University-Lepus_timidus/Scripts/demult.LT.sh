#!/bin/bash

#SBATCH --job-name=radtags
#SBATCH --account=project_2002674
#SBATCH --time=72:00:00
#SBATCH --mem-per-cpu=200M
#SBATCH --cpus-per-task=8
#SBATCH --partition=small
#SBATCH --mail-type=END

###Array setup here
##SBATCH --array=1-80%5
##SBATCH --open-mode=truncate

#SBATCH --output=proc_radtag_LT_output_%j.txt
#SBATCH --error=proc_radtag_LT_errors_%j.txt


#Process_radtags
export PATH="/projappl/project_2002674/stacks-2.68/bin:$PATH"

#Variables
vault=/scratch/project_2002674/zsofia/Hybrid_hares/00_Data/raw_fastq/
outdir=/scratch/project_2002674/zsofia/Hybrid_hares/01_tagproc/out/
barcodir=/scratch/project_2002674/zsofia/Hybrid_hares/00_Data/barcode/

while read id; do
    f1=$vault/${id}_1.fq.gz
    f2=$vault/${id}_2.fq.gz
    bcd=$(echo $id | cut -d'_' -f1-3)
    bcdf=$bcd.txt
    echo "Running ${f1} and ${f2} with ${bcdf} as barcode file"    

    mkdir -p $outdir/$id/
    
    process_radtags \
	-1 $f1 \
	-2 $f2 \
	-o $outdir/$id/ \
	-b $barcodir/$bcdf \
	-c \
	-q \
	-r \
	--renz_1 sbfI --renz_2 mspI \
	--adapter_1 GATCGGAAGAGCGGTTCAGCAGGAATGC --adapter_2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
	-E phred33 \
	--adapter_mm 1 \
	--threads 8 \
	--basename $id
	-w 0.1 -s 20 -t 143 \
	--len_limit 143 \
done < raw_id.LT.list #List of prefixes of multiplex fastqs
