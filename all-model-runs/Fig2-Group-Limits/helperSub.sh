#!/bin/bash
#SBATCH --job-name=redogrouplim
#SBATCH --account=fc_sirmodel
#SBATCH --partition=savio3
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=2
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=1
#SBATCH --time=72:00:00
module load r/3.6.3
module load r-packages
ht_helper.sh -m "r" -t RunScriptsSub.sh -p 2
