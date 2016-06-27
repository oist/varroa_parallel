#!/bin/bash

#SBATCH --job-name=switch_error
#SBATCH --partition=compute
#SBATCH --time=7-0
#SBATCH --mem-per-cpu=10G

### Calculate switch error for different filtering levels.
#   Phased data is in beagle/XX/ where XX is the filtering level.
#   Corresponding pre-phased drone data is in drones/fake_diploid.vcf.gz file.
#   Respective output is saved in drones/XX

# sbatch --array=1-16 --output=log/switch_error_X_%a.out --error=log/switch_error_X_%a.err switch_error.sh X
# where X is the filtering level 0..7

module load zlib/1.2.8

missing=( 99 95 90 85 80 70 60 50 )
ms=${missing[$1]}
chrom=$SLURM_ARRAY_TASK_ID

fin="beagle/"$ms"/Group"$chrom".vcf.gz"
base="Group"$chrom
fdrone="drones/fake_diploid.vcf.gz"

vcftools --gzdiff $fdrone --gzvcf $fin --keep drones/ids.txt --diff-switch-error --out "drones/""$ms"/""$base --chr LG$chrom

vcftools --gzvcf $f --keep drones/ids.txt --freq --out "drones/""$ms"/""$base --chr LG$chrom
