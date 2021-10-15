#!/bin/bash
#SBATCH --job-name=group-lim
#SBATCH --account=fc_sirmodel
#SBATCH --partition=savio3
#SBATCH --ntasks-per-node=6
#SBATCH --nodes=1
#SBATCH --time=48:00:00
module load r/4.0.3
module load r-packages
ht_helper.sh -m "r" -t RunScripts.sh -p 6
