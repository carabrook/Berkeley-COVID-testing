#!/bin/bash
#SBATCH --job-name=multigroup
#SBATCH --account=fc_sirmodel
#SBATCH --partition=savio3
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=7
#SBATCH --ntasks=7
#SBATCH --cpus-per-task=1
#SBATCH --time=48:00:00
module load r/3.6.3
module load r-packages
ht_helper.sh -m "r" -t RunScripts.sh -p 7