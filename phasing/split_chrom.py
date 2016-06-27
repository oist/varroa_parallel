import pdb
import sys
import os

"""
Generate vcfs of linkage groups in correct order.
../vcf/filtered/gatk_XX.vcf - input files, resulting from variant calling pipeline and filtering
For each filtering level XX a chrom/XX folder is created with a variant file of each linkage group in the correct ordering.
"""

missing=[99, 95, 90, 85, 80, 70, 60, 50]

# read group sizes
sizes = dict([i.rstrip().split() for i in open("sizes.txt")])

chromosomes = {}
for line in open("orientation.txt"):
	line=line.rstrip().split()
	chrom, group = line[0].split(".")
	chrom = chrom.replace("Group","")
	if chrom not in chromosomes:
		chromosomes[chrom] = []
	chromosomes[chrom].append((line[0],line[1]))

for ms in missing:
	os.mkdir("chrom/{}".format(ms))
	for chrom in chromosomes:
		for group, orientation in chromosomes[chrom]:
			# pdb.set_trace()
			if group == "Group2.12":
			# This group is not placed correctly according to Greg Hunt, so we'll rename it
				group2 = "Group1.19.5"
			else:
				group2 = group
			command = "vcftools --vcf ../vcf/filtered/gatk_{}.vcf --recode --chr {} -c | grep -v ^# ".format(ms, group)
			if orientation == "-":
			# reverse chromosome order and numbering of positions
				command +=  ("| tac | awk -v OFS=\"\\t\" -v size=%s \'{$2=size-$2; print}\'" %sizes[group])
			command += " | bgzip > chrom/{}/{}.vcf.gz; tabix -p vcf chrom/{}/{}.vcf.gz".format(ms,group2,ms,group2)
			print command
