#!/bin/bash

#SBATCH --time=10-00:00:00
#SBATCH --partition=regular
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=10000
#SBATCH --job-name=modelsele3
#SBATCH --output=modelsele3.log

module load Python/3.6.4-foss-2018a
srun python3 ModelSeleTP_DR_NH_clusterE1.py
