#!/bin/bash

#SBATCH --time=10-00:00:00
#SBATCH --partition=regular
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=30000
#SBATCH --job-name=MS2wfmTL
#SBATCH --output=MSBS2wfmTL.log

module load Python/3.6.4-foss-2018a
srun python3 BWMS2w_fm_TL.py
