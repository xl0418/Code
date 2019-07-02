#!/bin/bash

timescaling=(20000 40000 60000 80000 100000)
heri_sqr=(1 0.5)

for((t = 0; t <= 4; t++))
do
for((f = 0; f <= 1; f++))
do
echo "#!/bin/bash" > est_contrast$t$f
echo "#SBATCH --time=10-00:00:00" >> est_contrast$t$f
echo "#SBATCH --ntasks=1" >> est_contrast$t$f
echo "#SBATCH --nodes=1" >> est_contrast$t$f
echo "#SBATCH --cpus-per-task=20" >> est_contrast$t$f

echo "srun python3 BaleenWhaleEstContrastCluster.py ${timescaling[t]} ${heri_sqr[f]} " >> est_contrast$t$f
echo "rm est_contrast$t$f" >> est_contrast$t$f

sbatch --partition=gelifes --mem=60GB --job-name=BWcontrast$t$f --output=BWcontrast$t$f.log est_contrast$t$f

done
done