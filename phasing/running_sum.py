import sys,pdb

"""
Take GroupX.Y as input and give the length of all the scaffolds preceding it.
Requires file sizes.txt.
"""

# read group sizes
sizes = {}
for i in open("sizes.txt"):
	chrom, size = i.rstrip().split()
	chrom = chrom.split(".")
	scaf = ".".join(chrom[1:])
	chrom = chrom[0]
	if chrom not in sizes:
		sizes[chrom] = {}
	sizes[chrom][float(scaf)]=int(size)

chrom = sys.argv[1].split(".")
scaf = ".".join(chrom[1:])
chrom = chrom[0]

total = 0
keys = sorted(sizes[chrom].keys())
index = keys.index(float(scaf))
for i in keys[:index]:
	total += sizes[chrom][i]

print total
