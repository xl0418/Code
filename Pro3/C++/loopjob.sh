#!/bin/bash
phiarray=(0 0.001 0.01 0.1 1)
psiarray=(0 0.001 0.01 0.1 1)
for i in {0..4}
	do for j in {0..4}
		phivalue=${phiarray[i]}
		psivalue=${psiarray[j]}
		do
		 sbatch jc.sh $phivalue $psivalue $i $j
		done
	done