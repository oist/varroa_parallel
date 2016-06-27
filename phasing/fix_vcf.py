import sys
for line in open(sys.argv[1]):
	if line[0] == "#":
		print line,
		continue
	line = line.rstrip().split("\t")
	for i in range(9,len(line)):
		elements = line[i].split(":")
		if elements[0] == ".":
			line[i] = "./."
		else:
			for gl in elements[1].split(","):
				if gl == ".":
					line[i] = "./."
					break
	print "\t".join(line)
