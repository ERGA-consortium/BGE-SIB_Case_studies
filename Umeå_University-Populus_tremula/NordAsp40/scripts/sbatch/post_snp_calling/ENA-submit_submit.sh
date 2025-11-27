#!/bin/bash -l

#SBATCH -A hpc2n2025-120	
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 1
#SBATCH -t 2-00:00:00
#SBATCH -J ENA_submit
#SBATCH --output=/home/m/mimmie/projects/aspen_hpc2nstor2025-059/mimmi/NordAsp_40/reports/sbatch/ENA_submit/sbatch_R-%x_%j.out
#SBATCH --error=/home/m/mimmie/projects/aspen_hpc2nstor2025-059/mimmi/NordAsp_40/reports/sbatch/ENA_submit/sbatch_R-%x_%j.err

ml Java/17.0.6

cd /home/m/mimmie/projects/aspen_hpc2nstor2025-059/data/populus_tremula/fastq/NordAsp/raw/reads

ena_cli="/home/m/mimmie/software/ENA_webin-cli_v9.0.1/webin-cli-9.0.1.jar"

ena_user=$(cat ena_account.txt | cut -f1)
ena_password=$(cat ena_account.txt | cut -f2)

for manifest in $(ls /home/m/mimmie/projects/aspen_hpc2nstor2025-059/mimmi/NordAsp_40/ENA_submission/*.txt)
do 
    echo "Sumbitting "$(basename ${manifest})
    java -Xms2G -jar ${ena_cli} -context reads -username=${ena_user} -password=${ena_password} -manifest ${manifest} -sampleUpdate -submit
    echo "Done with "$(basename ${manifest})
    sleep 10m
done