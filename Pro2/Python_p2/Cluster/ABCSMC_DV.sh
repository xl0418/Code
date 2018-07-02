#!/bin/bash
#SBATCH --time=10-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --job-name=ABCSMC_DV
#SBATCH --output=ABCSMC_DV.log
#SBATCH --mem=5000
module load Tkinter/3.6.4-foss-2018a-Python-3.6.4
python ABC_SMC_run_tree.py
