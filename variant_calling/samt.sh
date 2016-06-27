#!/bin/bash

#SBATCH --job-name=samt
#SBATCH --partition=compute
#SBATCH --mem-per-cpu=10G
#SBATCH --output=log/samt_%a.out
#SBATCH --error=log/samt_%a.err

### Variant calling using samtools.
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

intfile="varroa/samt_"$SLURM_ARRAY_TASK_ID".bed"
rm -f $intfile
for s in ${rs[@]}; do
    ss=(${s//:/ })
    nm=${ss[0]}
    s=${ss[1]}
    ss=(${s//-/ })
    i1=${ss[0]}
    i2=${ss[1]}
    i1=$((i1 - 1))
    echo "$nm $i1 $i2" >> $intfile
done

fout="vcf/samt_"$SLURM_ARRAY_TASK_ID".vcf"

samtools mpileup -ugf $ref --positions $intfile --bam-list "bam_list" | bcftools call -vc -Ov > $fout

rm $intfile
