#!/bin/bash -l

#SBATCH -A hpc2n2025-120	
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 10
#SBATCH -t 3-00:00:00
#SBATCH -J stairwayplot2
#SBATCH --output=/home/m/mimmie/projects/aspen_hpc2nstor2025-059/mimmi/NordAsp_40/reports/sbatch/stairwayplot2/sbatch_R-%x_%j.out
#SBATCH --error=/home/m/mimmie/projects/aspen_hpc2nstor2025-059/mimmi/NordAsp_40/reports/sbatch/stairwayplot2/sbatch_R-%x_%j.err

ml Java/17.0.6

cd /proj/nobackup/hpc2nstor2025-059/mimmi/NordAsp_40
source .venv/bin/activate

# make stairwayplot2 input file i.e. blueprint file
python3 scripts/python/build_stairway_bluprint.py \
    -i output_data/vcf_filtered/NordAsp_biallelic_LDpruned.vcf.gz \
    -r "NordAsp40_r2" \
    -l 361795894 \
    -m 0.0000000025 \
    -g 15 \
    -d "output_data/stairway/" \
    -ds "/home/m/mimmie/software/stairway_plot_v2.2/" \
    -o "output_data/stairway/NordAsp_40_r2.blueprint" 

# configure stairwayplot2
java -cp /home/m/mimmie/software/stairway_plot_v2.2/stairway_plot_es Stairbuilder output_data/stairway/NordAsp_40_r2.blueprint

# run stairwayplot2
bash output_data/stairway/NordAsp_40_r2.blueprint.sh