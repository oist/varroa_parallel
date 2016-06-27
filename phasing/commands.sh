#!/bin/bash

#SBATCH --job-name=split
#SBATCH --partition=compute
#SBATCH --mem-per-cpu=4G
#SBATCH --time=2-0
#SBATCH --output=log/split_%a.out
#SBATCH --error=log/split_%a.err

old_IFS=$IFS
IFS=$'\n'
a=(`cat commands.txt`)
#a=(`cat coord_cmds.txt`)
IFS=$old_IFS
f=${a[$SLURM_ARRAY_TASK_ID]}

echo $f
eval $f

# sbatch --array=0-337 commands.sh
# sbatch --array=0-2703 commands.sh
