#!/bin/bash
#SBATCH --time=6-20:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --job-name=MCMCDV4wchain
#SBATCH --output=MCMCDV4wchain.log
#SBATCH --mem=1GB
module load Tkinter/3.6.4-foss-2018a-Python-3.6.4
python ABCMCMCchain4w.py
