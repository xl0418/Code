#!/bin/bash

phi=(0 0.001 0.01 0.1 1)
psi=(0 0.001 0.01 0.1 1)

for((s = 0; s <= 4; s++))
do
for((m = 0; m <= 4; m++))
do
echo "#!/bin/bash" > zMLjob$s$m
echo "#SBATCH --time=05:59:00" >> zMLjob$s$m
echo "./jc phi=${phi[s]} psi=${psi[m]} sA=5 sB=1 ticks=10e6 seed=1 file="test${phi[s]}${psi[m]}.m" " >> zMLjob$s$m
echo "rm zMLjob$s$m" >> zMLjob$s$m

#sbatch --partition=regular --mem=12GB --job-name=ML$s --mail-type=FAIL,TIME_LIMIT --mail-user=glaudanno@gmail.com zMLjob$s
sbatch --partition=regular --mem=1GB --job-name=JC$s$m --output=JC$s$m.log --mail-type=FAIL,TIME_LIMIT --mail-user=xl0418@gmail.com zMLjob$s$m

done
done
