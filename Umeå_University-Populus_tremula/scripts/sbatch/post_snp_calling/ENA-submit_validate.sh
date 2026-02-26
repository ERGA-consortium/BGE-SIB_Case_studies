#!/bin/bash -l

#SBATCH -A hpc2n	
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 1
#SBATCH -t 1-00:00:00
#SBATCH -J ENA_validate
#SBATCH --output=reports/sbatch/ENA_validate/sbatch_R-%x_%j.out
#SBATCH --error=reports/sbatch/ENA_validate/sbatch_R-%x_%j.err

ml Java/17.0.6

cd fastq/NordAsp/raw/reads/

ena_cli="software/ENA_webin-cli_v9.0.1/webin-cli-9.0.1.jar"

ena_user=$(cat ena_account.txt | cut -f1)
ena_password=$(cat ena_account.txt | cut -f2)

for manifest in $(ls ENA_submission/*.txt)
do 
    echo "Validating "$(basename ${manifest})
    java -Xms2G -jar ${ena_cli} -context reads -username=${ena_user} -password=${ena_password} -manifest ${manifest} -sampleUpdate -validate
    echo "Done with "$(basename ${manifest})
    sleep 5m
done

