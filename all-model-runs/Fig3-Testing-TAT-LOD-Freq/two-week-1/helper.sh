#!/bin/bash
#SBATCH --job-name=2wktest
#SBATCH --account=fc_sirmodel
#SBATCH --partition=savio3
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=9
#SBATCH --ntasks=9
#SBATCH --time=72:00:00
module load r/3.6.3
module load r-packages
ht_helper.sh -m "r" -t RunScripts.sh -p 9
