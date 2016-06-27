import os, pdb

"""
Generate a map from original coordinates to continuous coordinates.
Requires sizes.txt and orientation.txt files.
Coordinates are saved in coords folder.
"""

sizes = dict([i.rstrip().split() for i in open("sizes.txt")])

chromosomes = {}
for line in open("orientation.txt"):
	line=line.rstrip().split()
	chrom, group = line[0].split(".")
	chrom = chrom.replace("Group","")
	if chrom not in chromosomes:
		chromosomes[chrom] = []
	chromosomes[chrom].append((line[0],line[1]))

outfile = open(group+".vcf","w")

for chrom in chromosomes:
	for group, orientation in chromosomes[chrom]:
		# pdb.set_trace()
		if group == "Group2.12":
			# This group is not placed correctly according to Greg Hunt, so we'll rename it
			group2 = "Group1.19.5"
		else:
			group2 = group
		offset = os.popen("python running_sum.py %s" % group2).read().rstrip()
		# pdb.set_trace()
		command = "vcftools --vcf ../vcf/filtered/gatk.vcf --recode --chr {} -c | grep -v ^# ".format(group)
		if orientation == "-":
			command +=  ("| tac | awk -v OFS=\"\\t\" -v offset=%s -v size=%s \'{pos2=size-$2; print $1,$2,pos2+offset}\'" % (offset, sizes[group]))
		else:
			command +=  ("| awk -v offset=%s -v OFS=\"\\t\" \'{print $1,$2,$2+offset}\'" % offset)
		command += " > coords/{}.txt".format(group2)
		print command
