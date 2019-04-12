#!/bin/bash

timescaling=(20000 40000 80000)
dividing=(1 4)
fixm=(0 1)

for((t = 0; t <= 2; t++))
do
for((d = 0; d <= 1; d++))
do
for((f = 0; f <= 1; f++))
do
echo "#!/bin/bash" > zMLjob$t$d$f
echo "#SBATCH --time=10-00:00:00" >> zMLjob$t$d$f
echo "#SBATCH --ntasks=1" >> zMLjob$t$d$f
echo "#SBATCH --nodes=1" >> zMLjob$t$d$f
echo "#SBATCH --cpus-per-task=24" >> zMLjob$t$d$f

echo "srun python3 BWMS2.py "${timescaling[t]}" "${dividing[d]}" "${fixm[f]}" " >> zMLjob$t$d$f
echo "rm zMLjob$t$d$f" >> zMLjob$t$d$f

#sbatch --partition=regular --mem=12GB --job-name=ML$s --mail-type=FAIL,TIME_LIMIT --mail-user=glaudanno@gmail.com zMLjob$s
sbatch --partition=regular --mem=3GB --job-name=BWMS$t$d$f --output=BWMS$t$d$f.log zMLjob$t$d$f

done
done
done