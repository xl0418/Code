#!/bin/bash

#SBATCH --time=00:30:00
#SBATCH --partition=short
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=24
#SBATCH --mem-per-cpu=100
#SBATCH --job-name evopy

module load Python
srun python3 abcpp/abcpp/evo.py

