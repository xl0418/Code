#!/bin/bash
#SBATCH --time=12-20:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --job-name=MCMC_DV10w
#SBATCH --output=MCMC_DV10w.log
#SBATCH --mem=1GB
module load Tkinter/3.6.4-foss-2018a-Python-3.6.4
python ABC_MCMC_10w.py
