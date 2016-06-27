#!/bin/bash

#SBATCH --job-name=bays
#SBATCH --partition=largemem
#SBATCH --mem-per-cpu=12G
#SBATCH --time=2-0
#SBATCH --output=log/bays.out
#SBATCH --error=log/bays.err

### Filtering of varian calling results by freebayes, GATK and samtools
#   log - folder with output and error logs
#   vcf/filtered - output folder
### Results of GATK alignment after this filtering are used in furher analysis.

CALLERS=( freebayes gatk samt )

fout="vcf/filtered/all"

a=( "${CALLERS[@]/#/--vcf vcf/}" )
a=( "${a[@]/%/.vcf}" )
i=${a[@]}

baysic.pl --vcfOutFile $fout".vcf" --statsOutFile $fout".stats" --pvalCutoff 0.8 --countsOutFile $fout".cts" $i
