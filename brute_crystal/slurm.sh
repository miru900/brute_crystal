#!/bin/bash -i
# 
#SBATCH --job-name=IPCSP
#SBATCH --output=slurm.out
#SBATCH --nodelist=node8
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=36
#SBATCH --time=24:00:00
#SBATCH -p 36core_partition

python central.py
