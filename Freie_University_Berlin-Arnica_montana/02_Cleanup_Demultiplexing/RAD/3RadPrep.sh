#!/bin/bash
#SBATCH --job-name=3RadPrep           # Job name
#SBATCH --output=3RadPrep_output.txt  # Standard output file
#SBATCH --error=3RadPrep_error.txt    # Standard error file
#SBATCH --nodes=1                     # Number of nodes
#SBATCH --cpus-per-task=16            # Number of CPU cores per task
#SBATCH --mem=2G                      # Memory [K|M|G|T]
#SBATCH --time=3:00:00                # Maximum runtime (D-HH:MM:SS)
#SBATCH --qos=standard                # Quality of service
#SBATCH --mail-type=END               # Send email at job completion
#SBATCH --mail-user=reichel@zedat.fu-berlin.de    # Email address for notifications

# say hi
echo "Starting Job $SLURM_ARRAY_TASK_ID on $HOSTNAME"

# setup work environment
module load Miniconda3
eval "$(mamba shell hook --shell bash)"
cd /home/reichel/tools/rad/radseq-preprocessing-pipeline-main/
mamba activate radPreprocessing

# run job
snakemake --cores 16 --use-conda --conda-prefix scratch/reichel/ArnicaSNP/03_Preprocessed --show-failed-logs --rerun-incomplete
# further options
# --show-failed-logs --rerun-incomplete

#say bye
echo "Job $SLURM_ARRAY_TASK_ID on $HOSTNAME finished"