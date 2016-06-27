#!/bin/bash

#SBATCH --job-name=mrg
#SBATCH --partition=compute
#SBATCH --mem-per-cpu=4G
#SBATCH --time=2-0

### Pool the linkage groups into chromosomes. This script uses running_sum.py to re-number the positions continuously along an entire linkage group. It also uses fix_vcf.py to fix coding of missing genotypes, so that they are compatible with BEAGLE.
# For each filtering level XX a linkage groups in the folder chrom/XX folder are pooled into a respective folder merged/XX.


missing=( 99 95 90 85 80 70 60 50 )
ms=${missing[$1]}

chrom=$SLURM_ARRAY_TASK_ID
tmpfile=merged/"$ms"/Group"$chrom".vcf
outfile=merged/"$ms"/Group"$chrom".vcf.gz

head -5675 ../vcf/filtered/gatk_"$ms".vcf.gz > $tmpfile

for i in `ls -1 chrom/"$ms"/Group"$chrom".*.vcf.gz | sort -V` ; do
	echo $i
	# change chromosomes and positions to have continuous numbering of bases
	group=`echo $i | sed 's|^chrom/'$ms'/||' | sed 's/.vcf.gz//'`
	# count how many bases to add to make the chromosome continuously numbered
	echo $group
	total=`python running_sum.py $group`
	zcat $i | awk -v OFS="\t" -v chrom=$chrom -v total=$total '{$1="LG"chrom; $2+=total; print}' >> $tmpfile
done

python fix_vcf.py $tmpfile | bgzip > $outfile

rm $tmpfile

tabix -p vcf $outfile
