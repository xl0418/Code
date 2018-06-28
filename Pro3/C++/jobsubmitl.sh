#!/bin/bash
#SBATCH --time=05:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=1g
./jc phi=$1 psi=$2 sA=5 sB=1 ticks=10e6 seed=1 file="test$3.m"

