#!/bin/bash
#SBATCH --time=00:01:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --job-name=jcsim_%i_%j
#SBATCH --output=jcsim_%i_%j.log
#SBATCH --mem=100m
echo "phi value $d: $.4f\n" $i $phivalue

echo "psi value $d: $.4f" $j $psivalue

#./jc phi=$phivalue psi=$psivalue sA=5 sB=1 ticks=10e7 seed=1 file="test_phi$phivalue_psi.m"
