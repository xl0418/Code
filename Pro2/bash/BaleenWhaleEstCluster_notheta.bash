#!/bin/bash

timescaling=(20000 40000 60000 80000 100000)
heri_sqr=(1 0.5)

for((t = 0; t <= 4; t++))
do
for((f = 0; f <= 1; f++))
do
echo "#!/bin/bash" > est_notheta$t$f
echo "#SBATCH --time=5-00:00:00" >> est_notheta$t$f
echo "#SBATCH --ntasks=1" >> est_notheta$t$f
echo "#SBATCH --nodes=1" >> est_notheta$t$f
echo "#SBATCH --cpus-per-task=20" >> est_notheta$t$f

echo "srun python3 BaleenWhale_Est_cluster_notheta.py ${timescaling[t]} ${heri_sqr[f]} " >> est_notheta$t$f
echo "rm est_notheta$t$f" >> est_notheta$t$f

sbatch --partition=gelifes --mem=60GB --job-name=BWest1$t$f --output=BWest1$t$f.log est_notheta$t$f

done
done