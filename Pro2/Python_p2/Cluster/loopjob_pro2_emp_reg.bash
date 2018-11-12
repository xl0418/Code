#!/bin/bash
#SBATCH --partition=regular

#gamma=(0 0.001 0.01 0.1 0.5 1)
#a=(0 0.001 0.01 0.1 0.5 1)

for((s = 0; s <= 5; s++))
do
for((m = 0; m <= 5; m++))
do
echo "#!/bin/bash" > REGjobs$s$m
echo "#SBATCH --time=9-23:59:00" >> REGjobs$s$m
echo "#SBATCH --partition=regular" >> REGjobs$s$m
echo "#SBATCH --ntasks=1" >> REGjobs$s$m
echo "#SBATCH --nodes=1" >> REGjobs$s$m
echo "#SBATCH --cpus-per-task=16" >> REGjobs$s$m
echo "export OMP_NUM_THREADS="'$SLURM_CPUS_PER_TASK'  >> REGjobs$s$m
echo "module load Python/3.6.4-foss-2018a" >> REGjobs$s$m
echo "srun python3 evo_loop_empirical.py "$s" "$m"" >> REGjobs$s$m
echo "rm REGjobs$s$m" >> REGjobs$s$m

sbatch --partition=regular --mem-per-cpu=1GB --job-name=empismc$s$m --output=empiREGsmc$s$m.log --mail-type=FAIL,TIME_LIMIT --mail-user=xl0418@gmail.com REGjobs$s$m

done
done
