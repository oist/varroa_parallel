from glob import glob
import re,pdb

"""
fix tped for selscan
 - add map info
 - change values to zeros and ones
 - add chromosome number
"""

cMb = 19.0/10**6 # centimorgans per base

for f in glob("tped/*tped"):
	infile = open(f)
	outfile = open(f.replace(".",".fixed."),"w")
	group = re.search(r'[0-9]+', f).group(0) # extract group info
	for line in infile:
		line = line.rstrip().split()
		line[0] = group
		line[2] = str(int(line[1].split(":")[1])*cMb)  # estimate map distance
		ref = line[4]
		for i in range(4,len(line)):
			if line[i] == ref:
				line[i] = "0"
			else:
				line[i] = "1"
		outfile.write(" ".join(line)+"\n")
	outfile.close()
	infile.close()
