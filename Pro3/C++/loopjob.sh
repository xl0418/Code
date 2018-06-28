#SBATCH --array=1-25
INPUTFILE=parameterfile.txt
# get n-th line from $INPUTFILE
ARGS=$(cat $INPUTFILE | head -n $SLURM_ARRAY_TASK_ID | tail -n 1)
jobsubmitl.sh $ARGS