#!/bin/bash

#SBATCH --job-name=mrg_flt
#SBATCH --partition=compute
#SBATCH --mem-per-cpu=8G
#SBATCH --time=2-0
#SBATCH --output=log/flt_%a.out
#SBATCH --error=log/flt_%a.err

### Merging of previously split drone and worker data after missing data-based filtering on variants in workers.
#   log - folder with output and error logs
#   vcf/filtered/gatk_workers.vcf, vcf/filtered/gatk_fakedrones.vcf - previously created variant file based on BAYSIC results
#   vcf/filtered/gatk_workers_XX.kept.sites - output files of the filtering based on different levels of missing data
#   vcf/filtered/gatk_XX.vcf.gz - output files

missing=( 99 95 90 85 80 70 60 50 )
ms=${missing[$SLURM_ARRAY_TASK_ID]}

f_sites="vcf/filtered/gatk_workers_"$ms".kept.sites"

f_drones="vcf/filtered/gatk_fakedrones.vcf"
f_workers="vcf/filtered/gatk_workers.vcf"

fout1="vcf/filtered/gatk_workers_"$ms"_tmp"
fout2="vcf/filtered/gatk_fakedrones_"$ms"_tmp"
fout="vcf/filtered/gatk_"$ms".vcf.gz"

vcftools --vcf $f_workers --recode --recode-INFO-all --positions $f_sites --out $fout1
vcftools --vcf $f_drones  --recode --recode-INFO-all --positions $f_sites --out $fout2

bgzip $fout1".recode.vcf"
bgzip $fout2".recode.vcf"

tabix -p vcf $fout1".recode.vcf.gz"
tabix -p vcf $fout2".recode.vcf.gz"

vcf-merge $fout1".recode.vcf.gz" $fout2".recode.vcf.gz" > $fout

rm $fout1".recode.vcf.gz"
rm $fout2".recode.vcf.gz"
