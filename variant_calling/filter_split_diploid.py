import sys

### From the results of variant calling and BAYSIC filtering variants where drone individuals are diploid are removed. Variants are split into workers and drone data for further analysis
#   log - folder with output and error logs
#   vcf/filtered - output folder
#   usage:
#	srun --job-name=filter --partition=compute --mem-per-cpu=8G --output=log/filter.out --error=log/filter.err python filter_split_diploid.py

def get_drone_list(fl="bam_list"):
	res = []
	f = open(fl, 'r')
	for line in f:
		t = line.split("/")
		if ("drones" in t):
			res.append(t[len(t)-1].split(".")[0])
	f.close()
	return res

drone_names=get_drone_list()
vcfs = 0
kept = 0
drone_is = []
other_is = []
drone_pairs = []
fin="vcf/filtered/gatk.vcf"
fout_workers="vcf/filtered/gatk_workers.vcf"
fw = open(fout_workers, 'w')
fout_drones="vcf/filtered/gatk_fakedrones.vcf"
fd = open(fout_drones, 'w')

for line in open(fin):
	line = line.rstrip()
	if (line.startswith("#CHROM")):
		t = line.split("\t")
		for drone_name in drone_names:
			drone_is.append(t.index(drone_name))
		sys.stderr.write(" ".join(map(str, drone_is)))
		other_is = []
		for i in range(9,len(t)):
			if (not (i in drone_is)):
				other_is.append(i)
		for i in range(0, len(drone_is)-1, 2):
			drone_pairs.append((drone_is[i], drone_is[i+1]))
		line1 = t[:9]
		for i in other_is:
			line1.append(t[i])
		line2 = t[:9]
		for (d1, d2) in drone_pairs:
			line2.append(t[d1]+"_"+t[d2])
		fw.write("\t".join(line1) + "\n")
		fd.write("\t".join(line2) + "\n")
	elif (line.startswith("#")):
		fw.write(line + "\n")
		fd.write(line + "\n")
	else:
		t = line.split("\t")
		haploid = True
		i = 0
		while (haploid and (i < len(drone_is))):
			s = t[drone_is[i]]
			if (s != "."):
				s = s.split(":")[0]
				if (s != "."):
					ss = s.split("/")
					s1 = ss[0]
					s2 = ss[1]
					if (s1 != s2):
						haploid=False
			i += 1
		if (haploid):
			line1 = t[:9]
			for i in other_is:
				line1.append(t[i])
			fw.write("\t".join(line1) + "\n")
			line2 = t[:9]
			for (d1, d2) in drone_pairs:
				el1 = t[d1].split(":")
				gt1 = el1[0].split("/")
				el2 = t[d2].split(":")
				gt2 = el2[0].split("/")
				s = ""
				if gt1[0] == "." or gt2[0] == ".":
					gt1 = [".","."]
				else:
					gt1[1] = gt2[0]
					s = ":" + ":".join(el1[1:])
				line2.append("/".join(gt1) + s)
			fd.write("\t".join(line2) + "\n")
			kept += 1
		vcfs += 1

fw.close()
fd.close()

sys.stderr.write("kept {0} out of {1}\n".format(kept, vcfs))
