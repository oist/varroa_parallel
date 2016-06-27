#!/bin/bash

#SBATCH --job-name=fbayes
#SBATCH --partition=compute
#SBATCH --mem-per-cpu=12G
#SBATCH --time=2-0
#SBATCH --output=log/fb_%a.out
#SBATCH --error=log/fb_%a.err

### Variant calling using freebayes.
#   log - folder with output and error logs
#   vcf - output alignment folder

ref="ref/Amel_4.5.AGP.linearScaffold.fa"

getArray() {
    local i=0
    array=()
    while IFS= read -r line
    do
        array+=( "$line" )
    done < "$1"
}
getArray "ref/refs.txt"
rs=( ${array[$SLURM_ARRAY_TASK_ID]} )

fout="vcf/freebayes_"$SLURM_ARRAY_TASK_ID".vcf"
rm -f $fout
for r in ${rs[@]}
do
    freebayes --theta 0.002 -f $ref -r $r --bam-list "bam_list" --use-best-n-alleles 6 >> $fout
done
