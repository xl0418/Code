#!/bin/bash
#SBATCH --time=5-20:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --job-name=ABCMCMC
#SBATCH --output=ABCMCMC.log
#SBATCH --mem=5000
module load Tkinter/3.6.4-foss-2018a-Python-3.6.4
python ABC_MCMC_run_cluster.py
