#!/bin/bash

#SBATCH --job-name=assoc
#SBATCH --partition=compute
#SBATCH --time=7-0
#SBATCH --mem-per-cpu=90G

### Run phasing for the reoriented chromosome variants.
# For each filtering level XX phased chromosome data is created in beagle/XX.
# sbatch --array=1-16 --output=log/assoc_X_%a.out --error=log/assoc_X_%a.err assoc.sh X
# where X is the filtering level 0..7

missing=( 99 95 90 85 80 70 60 50 )
ms=${missing[$1]}

chrom=$SLURM_ARRAY_TASK_ID

fin=merged/"$ms"/Group"$chrom".vcf.gz
fout=beagle/"$ms"/Group"$chrom"

java -Xmx90g -jar beagle.14Jan16.841.jar gt=$fin out=$fout
