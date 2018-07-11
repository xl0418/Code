#!/bin/bash
#SBATCH --time=5-20:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --job-name=Cali_1
#SBATCH --output=Cali_1.log
#SBATCH --mem=1GB
module load Tkinter/3.6.4-foss-2018a-Python-3.6.4
python Calibration1_cluster.py
