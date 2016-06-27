import  gzip
import pdb

"""
Merge pre-phased drone variants into synthetic diploids that will be used in further quality checks. Reorder the chromosome positions.
../vcf/filtered/gatk.vcf - variants resulting from variant calling pipeline and filtering

python merge_drones.py | bgzip > drones/fake_diploid.vcf.gz; tabix -p vcf drones/fake_diploid.vcf.gz
"""

drones = [ "159", "160", "64", "65", "66", "67", "68", "69", "70", "71", "72", "73", "74", "75", "D1", "D2", "D3", "D4", "D5" ]
drone_ps = []

# read mappings from old to new coordinates
coords = {}
for line in open("coords/all.txt"):
	group, orig, new = line.rstrip().split()
	lg = group.split(".")[0].replace("Group","LG")
	coords[(group,orig)] = (lg,new)

for line in open("../vcf/filtered/gatk.vcf"):
	line = line.rstrip();
	if line[0] == "#":
		if line.startswith("#CHROM"):
			line = line.split()
			outstr = "\t".join(line[:9])
			prev_i = -1
			for i in range(9,len(line)):
				if (line[i] in drones):
					if prev_i == -1:
						prev_i = i
					else:
						drone_ps.append((prev_i, i))
						outstr += "\t"+line[prev_i]+"_"+line[i]
						prev_i = -1
			line = outstr
		print line
	else:
		line = line.split("\t")
	#switch drone coordinates to new mapping
        if (line[0],line[1]) in coords:
			line[0],line[1] = coords[(line[0],line[1])]
			line[8] = "GT"
			outstr = "\t".join(line[:9])

			for (i1, i2) in drone_ps:
				el1 = line[i1].split(":")
				gt1 = el1[0].split("/")
				el2 = line[i2].split(":")
				gt2 = el2[0].split("/")

				if gt1[0] == "." or gt2[0] == ".":
					gt1 = [".","."]
				else:
					gt1[1] = gt2[0]
				outstr += "\t" + "|".join(gt1)

			print outstr
