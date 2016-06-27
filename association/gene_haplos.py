import sys,gzip,os,operator,copy

"""
Extract individuals' haplotypes of genes based on variants. Save into fasta files for further fst analysis.

python gene_haplos.py <chromosome>
"""

chr = int(sys.argv[1])
pops = ["az", "tx"]
stgs = ["early", "late"]
grs = ["az_early", "az_late", "tx_early", "tx_late"]

fin="../phasing/beagle/90/Group%d.vcf.gz" % (chr)
fout_fasta="genes/fasta/haplotypes_%s_%s_%d.fasta"

def read_samples():
	res = {}
	for gr in grs:
		f = open(gr + ".txt", 'r')
		v = f.read()
		v = v.strip().split('\n')
		res[gr] = v
		f.close()
	return res

def read_genes():
	res = []
	f = open("genome/amel_reoriented.gff", 'r')
	for line in f:
		line = line.rstrip()
		t = line.split("\t")
		if ((t[0] == ("Group%d" % (chr))) and (t[2] == "gene")):
			nm = t[8]
			nm = nm[3:]
			res.append((nm, int(t[3]), int(t[4])))
	f.close()
	res = sorted(res, key=operator.itemgetter(1))
	print("genes: %d" % (len(res)))
	return res

def which_genes(genes, gene_i, pos):
	all_incl = False
	i = gene_i
	res = []
	while ((i < len(genes)) and (not all_incl)):
		if (pos >= genes[i][1]) and (pos <= genes[i][2]):
			res.append(genes[i][0])
			i += 1
		elif (pos > genes[i][2]):
			i += 1
		else:
			all_incl = True
	if (len(res) == 0):
		while ((gene_i < len(genes)) and (pos > genes[gene_i][2])):
			gene_i += 1
	return (res, gene_i)

def pos_refalt(t, inds):
	refalt = [t[3]]
	refalt.extend(t[4].split(','))
	maxlen = 0
	for s in refalt:
		maxlen = max(maxlen, len(s))
	for i in range(len(refalt)):
		refalt[i] = ("-" * (maxlen-len(refalt[i]))) + refalt[i]

	res = []
	for i in inds:
		s = t[i]
		s = s[0:3].split('|')
		res.append(refalt[int(s[0])])
		res.append(refalt[int(s[1])])
	return res

def update_genes(curr_genes, genes, v, pos):
	k = curr_genes.keys()
	k.extend(genes)
	all_genes = list(set(k))
	finished_haplo = {}
	v["pos"] = [pos]
	for gene in all_genes:
		if (gene in genes) and (not (gene in curr_genes.keys())):
			# new gene
			curr_genes[gene] = copy.deepcopy(v)
		elif (not (gene in genes)) and (gene in curr_genes.keys()):
			# gene haplotype finished
			finished_haplo[gene] = curr_genes[gene]
			del curr_genes[gene]
		else:
			# extend gene haplotype
			curr_v = curr_genes[gene]
			curr_v["pos"].append(pos)
			for gr in grs:
				gr_v = curr_v[gr]
				for i in range(len(gr_v)):
					gr_v[i] = gr_v[i] + v[gr][i]
				curr_v[gr] = gr_v
			curr_genes[gene] = curr_v
	return (curr_genes, finished_haplo)

def add_fasta(haplos):
	for gene in haplos.keys():
		gene_h = haplos[gene]
		for gr in grs:
			gr_haplo = gene_h[gr]
			gr_nms = gr_smpls[gr]
			for i in range(len(gr_nms)):
				h1 = gr_haplo[i*2]
				h1 = h1.replace('-', '')
				h2 = gr_haplo[i*2+1]
				h2 = h2.replace('-', '')
				fw_fasta[gr].write(">" + gene + "|" + gr_nms[i] + "|1\n")
				fw_fasta[gr].write(h1 +"\n")
				fw_fasta[gr].write(">" + gene + "|" + gr_nms[i] + "|2\n")
				fw_fasta[gr].write(h2 +"\n")

def remove_write_gene_fasta(haplos, remov):
	for gene in haplos.keys():
		gene_h = haplos[gene]
		pos = gene_h["pos"]
		dr = "fasta"
		if (remov):
			dr = "filtered_fasta"
		os.mkdir("genes/%s/%s" % (dr,gene))
		fw = open("genes/%s/%s/%s.fasta" % (dr,gene,gene), 'w')

		for pop in pops:
			gr_early = pop + "_early"
			hs_early = gene_h[gr_early]
			gr_late = pop + "_late"
			hs_late = gene_h[gr_late]

			gr_nms = gr_smpls[gr_early]
			for i in range(len(gr_nms)):
				h1 = hs_early[i*2]
				h2 = hs_early[i*2+1]
				fw.write(">" + gr_nms[i] + "|1\n")
				fw.write(h1 + "\n")
				fw.write(">" + gr_nms[i] + "|2\n")
				fw.write(h2 + "\n")

			gr_nms = gr_smpls[gr_late]
			for i in range(len(gr_nms)):
				h1 = hs_late[2*i]
				h2 = hs_late[2*i+1]
				if (remov):
					if (h1 in hs_early):
						fw.write(">" + gr_nms[i] + "|1\n")
						fw.write(h1 + "\n")
					if (h2 in hs_early):
						fw.write(">" + gr_nms[i] + "|2\n")
						fw.write(h2 + "\n")
				else:
					fw.write(">" + gr_nms[i] + "|1\n")
					fw.write(h1 + "\n")
					fw.write(">" + gr_nms[i] + "|2\n")
					fw.write(h2 + "\n")

		fw.close()

def remove_haplos_from_lines(haplos, pos_lines):
	# which are the novel haplotypes - remove them from the pos_lines
	for gene in haplos.keys():
		gene_h = haplos[gene]
		pos = gene_h["pos"]
		for pop in pops:
			gr_early = pop + "_early"
			hs_early = gene_h[gr_early]
			gr_late = pop + "_late"
			hs_late = gene_h[gr_late]
			for i in range(len(gr_inds[gr_late])):
				h1 = hs_late[2*i] in hs_early
				h2 = hs_late[2*i+1] in hs_early
				if h1 and h2:
					# novel haplotype
					ind = gr_inds[gr_late][i]
					for p in pos:
						pos_lines[p][ind] = ".|.:."
	return pos_lines

gr_smpls = read_samples()
genes = read_genes()
gr_inds = {}
gene_i = 0
finished_genes = 0
curr_genes = {}
pos_lines = {}
pos_genes = {}
for line in gzip.open(fin):
	line = line.rstrip()
	if (line.startswith("#CHROM")):
		t = line.split("\t")
		for gr in grs:
			nms = gr_smpls[gr]
			v = []
			for nm in nms:
				v.append(t.index(nm))
			gr_inds[gr] = v
	elif (not line.startswith("#")):
		t = line.split("\t")
		pos = int(t[1])
		(genes_in_pos, gene_i) = which_genes(genes, gene_i, pos)
		v = {}
		if (len(genes_in_pos) > 0):
			pos_lines[pos] = t
			pos_genes[pos] = genes_in_pos
			for gr in grs:
				v[gr] = pos_refalt(t, gr_inds[gr])
		(curr_genes, finished_haplo) = update_genes(curr_genes, genes_in_pos, v, pos)
		if (len(finished_haplo) > 0):
			finished_genes += len(finished_haplo)
			print("finished: %d gene_i: %d" % (finished_genes, gene_i))
			#add_fasta(finished_haplo)
			remove_write_gene_fasta(finished_haplo, False)
			### pos_lines = remove_haplos_from_lines(finished_haplo, pos_lines)
			for gene in finished_haplo.keys():
				# remove finished genes from pos_genes
				for p in finished_haplo[gene]["pos"]:
					pos_genes[p].remove(gene)
			for p in pos_genes.keys():
				if (len(pos_genes[p]) == 0):
					del pos_lines[p]
					del pos_genes[p]
