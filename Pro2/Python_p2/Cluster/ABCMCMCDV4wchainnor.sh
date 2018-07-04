#!/bin/bash
#SBATCH --time=7-20:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --job-name=MC4wnor
#SBATCH --output=MC4wnor.log
#SBATCH --mem=1GB
module load Tkinter/3.6.4-foss-2018a-Python-3.6.4
python ABCMCMC4wchainnor.py
