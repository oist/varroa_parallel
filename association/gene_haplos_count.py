import sys,gzip,os,operator,copy

"""
Count all genes' haplotypes in the AZ and TX pre- and pos-varroa populations. The counts are saved in csv tables with a row per haplotype and four columns corresponding to the populations.

python gene_haplos_count.py <chromosome>
"""

chr = int(sys.argv[1])
pops = ["az", "tx"]
stgs = ["early", "late"]
grs = ["az_early", "az_late", "tx_early", "tx_late"]

fin="phasing/Group%d.vcf.gz" % (chr)


def read_samples():
	res = {}
	for gr in grs:
		f = open("pop/" + gr + ".txt", 'r')
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


def get_all_haplos(gene_h):
    res = []
    for pop in ["az_early", "az_late", "tx_early", "tx_late"]:
        res.extend(gene_h[pop])
    return list(set(res))


def count_haplos(haplos):
    fw_out = open("haplos/count/chr%i.txt" % chr, 'a')
    for gene in haplos.keys():
		gene_h = haplos[gene]
		all_hs = get_all_haplos(gene_h)
		#pos = gene_h["pos"]
		for i in range(len(all_hs)):
			h = all_hs[i]
			az_e = gene_h["az_early"].count(h)
			az_l = gene_h["az_late"].count(h)
			tx_e = gene_h["tx_early"].count(h)
			tx_l = gene_h["tx_late"].count(h)
			fw_out.write("%s_%i,%i,%i,%i,%i\n" % (gene,i,az_e,az_l,tx_e,tx_l))
    fw_out.close()


gr_smpls = read_samples()
genes = read_genes()
gr_inds = {}
gene_i = 0
finished_genes = 0
curr_genes = {}
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
			for gr in grs:
				v[gr] = pos_refalt(t, gr_inds[gr])
		(curr_genes, finished_haplo) = update_genes(curr_genes, genes_in_pos, v, pos)
		if (len(finished_haplo) > 0):
			finished_genes += len(finished_haplo)
			print("finished: %d gene_i: %d" % (finished_genes, gene_i))
			count_haplos(finished_haplo)
