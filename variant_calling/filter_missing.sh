#!/bin/bash

#SBATCH --job-name=flt
#SBATCH --partition=compute
#SBATCH --mem-per-cpu=8G
#SBATCH --time=2-0
#SBATCH --output=log/flt_%a.out
#SBATCH --error=log/flt_%a.err

### Filtering out results of variant calling in workers based on the amount of missing data, allowing for 1, 5, 10, 15, 20, 30, 40 or 50% of missing individual data for every variant.
#   log - folder with output and error logs
#   vcf/filtered/gatk_workers.vcf - previously created variant file based on BAYSIC results
#   vcf/filtered - output folder

missing=( 99 95 90 85 80 70 60 50 )
ms=${missing[$SLURM_ARRAY_TASK_ID]}

fin="vcf/filtered/gatk_workers.vcf"
fout="vcf/filtered/gatk_workers_"$ms

vcftools --vcf $fin --kept-sites --minQ 20 --max-missing 0.$ms --out $fout
