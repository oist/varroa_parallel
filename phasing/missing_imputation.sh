#!/bin/bash

#SBATCH --job-name=miss
#SBATCH --partition=compute
#SBATCH --mem-per-cpu=90G
#SBATCH --time=5-0

### Calculate proportion of the missing variants of the not imputated in the imputation process.
#   Imputed data is in beagle/XX/ where XX is the filtering level.
#   Corresponding pre-phased drone data is in drones/fake_diploid.vcf.gz file.
#   Respective output is saved in imputation/missing_imputed_chrXX_Y.txt where XX is the chromosome and Y filtering level.

# sbatch --array=1-16 --output=log/missing_imp_%a_X.out --error=log/missing__imp_%a_X.err missing_imputation.sh X
# where X is the filtering level 0..7

python missing_imputation.py $SLURM_ARRAY_TASK_ID $1
