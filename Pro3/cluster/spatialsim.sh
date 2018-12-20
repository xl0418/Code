#!/bin/bash
#SBATCH --time=9-23:59:00
#sbatch --partition=regular
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=18
#sbatch --mem=12GB 
#sbatch --job-name=pro3sim
#sbatch --mail-type=FAIL,TIME_LIMIT 
#sbatch --mail-user=xl0418@gmail.com 
./jc batch=spatialpara.txt
