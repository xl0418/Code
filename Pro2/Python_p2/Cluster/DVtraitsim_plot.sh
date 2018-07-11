#!/bin/bash
#SBATCH --time=6-20:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --job-name=simplot
#SBATCH --output=simplot.log
#SBATCH --mem=1GB
module load Tkinter/3.6.4-foss-2018a-Python-3.6.4
python DVtraitsim_plot.py