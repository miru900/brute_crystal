#!/bin/bash -i
# 
#SBATCH --job-name=IPCSP
#SBATCH --output=slurm.out
#SBATCH --nodelist=node10
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --time=24:00:00
#SBATCH -p 32core_partition

python central.py
