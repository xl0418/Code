#!/bin/bash
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --job-name=jc
#SBATCH --output=jc.log
#SBATCH --mem=5000
jc phi=0.1 psi=0.1 sA=1000 sB=1000 ticks=10^6 seed=96 file="cluster.m"

