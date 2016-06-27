import gzip
import numpy
import sys
import pickle

drones = [ "159", "160", "64", "65", "66", "67", "68", "69", "70", "71", "72", "73", "74", "75", "D1", "D2", "D3", "D4", "D5" ]
drone_ps = [ "159_160", "64_65", "66_67", "68_69", "70_71", "72_73", "74_75", "D1_D2", "D3_D4"]

# read mappings from old to new coordinates
def read_coords():
    coords = {}
    for line in open("coords/all.txt"):
        group, orig, new = line.rstrip().split()
        lg = group.split(".")[0].replace("Group","LG")
        coords[(group,orig)] = (lg,new)
    return coords

def read_renumber(fname, samples, coords):
    res = {}
    for gr in range(1,17):
        res["LG%d" % (gr)] = {}
    sample_is = {}
    for line in open(fname):
        line = line.rstrip();
        if line[0] == "#":
            if line.startswith("#CHROM"):
                line = line.split()
                for i in range(9,len(line)):
                    if (line[i] in samples):
                        sample_is[line[i]] = i
        else:
            line = line.split("\t")
	#switch drone coordinates to new mapping
            if (line[0],line[1]) in coords:
                gr, c = coords[(line[0],line[1])]
                d = {}
                for samp in sample_is.keys():
                    i = sample_is[samp]
                    el = line[i].split(":")
                    d[samp] = el[0].split("/")[0]
                res[gr][c] = d
    return res

def drone2o():
    d2o = {}
    for s in drone_ps:
        ds = s.split("_")
        d2o[ds[0]] = ds[1]
        d2o[ds[1]] = ds[0]
    return d2o

def got_removed(ds, d2o):
    res = []
    for d in ds.keys():
        if (d != "D5") and (ds[d] != ".") and ds[d2o[d]] == ".":
            res.append(d)
    return res

def read_compare(chr, ms, prechr):
    fname = "beagle/%d/Group%d.vcf.gz" % (ms, chr)
    #prechr = premerge["LG%d" % chr]
    d2o = drone2o()
    sum_all = 0
    sum_wrong = 0
    sample_is = {}
    for line in gzip.open(fname):
        line = line.rstrip();
        if line[0] == "#":
            if line.startswith("#CHROM"):
                line = line.split()
                for i in range(9,len(line)):
                    if (line[i] in drone_ps):
                        sample_is[line[i]] = i
        else:
            line = line.split("\t")
            if not (line[1] in prechr.keys()):
                print("Missing site!! chr %d site %s" % (chr, line[1]))
            else:
                site = prechr[line[1]]
                rmd = got_removed(site, d2o)
                if len(rmd) > 0:
                    imputed = {}
                    for dp in sample_is.keys():
                        ds = dp.split("_")
                        i = sample_is[dp]
                        el = line[i].split(":")
                        t = el[0].split("|")
                        imputed[ds[0]] = t[0]
                        imputed[ds[1]] = t[1]
                    for d in rmd:
                        sum_all += 1;
                        if (site[d] != imputed[d]):
                            sum_wrong += 1
    return (sum_wrong, sum_all)

coords = read_coords()
premerge = read_renumber("../vcf/filtered/gatk.vcf", drones, coords)
 for chr in range(1,17):
    with open("imputation/premerge_chr%d.pickle" % chr, 'wb') as f:
        pickle.dump(premerge["LG%d" % chr], f, pickle.HIGHEST_PROTOCOL)

chr = int(sys.argv[1])
ms_i = int(sys.argv[2])

premerge = {}
with open("imputation/premerge_chr%d.pickle" % chr, 'rb') as f:
    premerge = pickle.load(f)

mss = [ 50, 60, 70, 80, 85, 90, 95, 99 ]
ms = mss[ms_i]
wrng, sum_all = read_compare(chr, ms, premerge)
l = "%d,%d,%d,%d,%f" % (ms, chr, wrng, sum_all, float(wrng)/sum_all)
with open("imputation/missing_imputed_chr%d_%d.txt" % (chr,ms), 'w') as f:
    f.write(l)

# srun --partition=compute --mem-per-cpu=90G -n 1 --output=log/miss.out --error=log/miss.err python missing_imputation.py &
