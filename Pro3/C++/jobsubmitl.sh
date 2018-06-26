#!/bin/bash
#SBATCH --time=05:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --job-name=jcsim
#SBATCH --output=jcsim.log
#SBATCH --mem=1g
./jc phi=0.1 psi=0.1 sA=5 sB=1 ticks=10e7 seed=1 file="test1.m"
./jc phi=0.1 psi=0.1 sA=5 sB=1 ticks=10e7 seed=2 file="test2.m"
./jc phi=0.1 psi=0.1 sA=5 sB=1 ticks=10e7 seed=3 file="test3.m"
