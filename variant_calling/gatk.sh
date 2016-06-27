#!/bin/bash

#SBATCH --job-name=gatk
#SBATCH --partition=compute
#SBATCH --mem-per-cpu=20G
#SBATCH --time=2-0
#SBATCH --output=log/gatk_%a.out
#SBATCH --error=log/gatk_%a.err

### Variant calling using GATK.
#   log - folder with output and error logs
#   vcf - output alignment folder

ref="ref/Amel_4.5.AGP.linearScaffold.fa"

a=(`ls data/alignments/workers/*.bam data/alignments/drones/*.bam`)
a=( "${a[@]/#/-I }" )
i=${a[@]}

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
intfile="varroa/gatk_"$SLURM_ARRAY_TASK_ID".intervals"
( IFS=$'\n'; echo "${rs[*]}" ) > $intfile
( IFS=$'\n'; echo "${rs[*]}" )

fout="vcf/gatk_"$SLURM_ARRAY_TASK_ID".vcf"

java  -Xmx30g -jar $GATK -nct 4 -allowPotentiallyMisencodedQuals -T HaplotypeCaller -R $ref -L $intfile \
        -hets 0.002 --max_alternate_alleles 6 --min_base_quality_score 20 \
        $i -o $fout

rm $intfile
