#!/bin/bash
#SBATCH --time=10-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --job-name=MS_ABCSMC
#SBATCH --output=MS_ABCSMC.log
#SBATCH --mem=5000
module load Tkinter/3.6.4-foss-2018a-Python-3.6.4
python ABCSMC_MS.py
