#!/usr/bin/env python

"""
Run evolutionary simulations of bee populations.


"""

import simuPOP as sim
from simuPOP import utils
import sys, types, os, math
import numpy as np
import random
import argparse
import pickle
from functools import partial

parser = argparse.ArgumentParser()
parser.add_argument('--corfile', type=str, help='output file for parameter values and corresponding correlation ceofficients')
parser.add_argument('--ms_file', type=str, help='temporary output file')
parser.add_argument('--selection', action='store_true', help='sample fitness parameter, if false only migration parameters will be used')
args = parser.parse_args()

# initial population size and relative size of the african population
# afr_size = 2.0 in all our simulations
global psize, afr_size

## gene bins is a pre-generated table with gene haplotype length and corresponding minor allele proportions observed in our data
## the values are binned according to the length of the haplotype into 20 bins
with open("gene_bins.pickle", "rb") as f:
    gene_bins = pickle.load(f)

## sampling 100 genes for each simulation
def sample_genes():
    res = []
    for i in range(len(gene_bins)):
        res.extend(random.sample(gene_bins[i], 5))
    return res

## function returning populations sizes for each generation
def get_sizes(gen):
    all_gen = [[0.659751037, 0.866666667],
               [0.651452282, 0.882352941],
               [0.543568465, 0.9],
               [0.31120332, 0.887640449],
               [0.049792531, 0.359550562],
               [0.244813278, 0.382022472],
               [0.406639004, 0.438202247],
               [0.35, 0.635416667],
               [0.336099585, 0.702702703],
               [0.360995851,0.544642857]]
    res = all_gen[gen-1]
    res = [ math.ceil(res[0]*psize), math.ceil(res[1]*psize), psize*afr_size ]
    return res

## generator returning individual's sex, 10% of individuals are male
def sex_func():
    i = 0
    while True:
        i = (i + 1) % 10
        if i >= 1:
            yield sim.MALE
        else:
            yield sim.FEMALE

## export the resutls of simulation into MS file format
## exported individuals are sampled in numbers equal to the numbers in our dataset
def export(pop, fl, append):
    subPops = ["AZ", "TX"]
    n_inds = { "AZ": (43, 87), "TX": (98, 102) }
    sim.stat(pop, popSize=True, subPops=["AZ", "TX"], numOfMales=True)
    n_fm = pop.dvars().numOfFemales

    count_pop = dict([(pop,0) for pop in subPops])

    with open(fl, 'a') as f:
        seg_sites = list(range(pop.numLoci(0)))
        if not append:
            f.write('\n//\nsegsites: %d\n' % len(seg_sites))
            f.write('positions: %s\n' % ' '.join([str(pop.locusPos(x)) for x in seg_sites]))

        for vsp in subPops:
            female_inds = []
            for ind in pop.individuals(vsp):
                if ind.sex() == sim.FEMALE:
                    female_inds.append(ind)
            for ind in random.sample(female_inds, n_inds[vsp][append]):
                count_pop[vsp] += 1
                for p in range(1): #2
                    geno = ind.genotype(p, 0)
                    f.write(''.join([str(geno[x]) for x in seg_sites]) + '\n')
    return count_pop

## simulation function
def simu(w, m1, m2, psize, afr_size, fl):
    print(fl)
    if os.path.exists(fl):
        os.remove(fl)

    matingScheme = sim.HaplodiploidMating(sexMode=sex_func, subPopSize=get_sizes)
    fitness = { (0,):1.0, (1,): 1.0, (0,0):w, (0,1):(1+w)/2, (1,1):1.0 }
    migrator = sim.BackwardMigrator(rate=[[0, m1, m2],
                                          [m1, 0, m2],
                                          [0, 0, 0]], begin=4, step=1)

    count_pop = {}
    sg = sample_genes()
    for i in range(len(sg)):
        (n_loci, gene_prop_az, gene_prop_tx) = sg[i]

        selector = sim.MlSelector([sim.MapSelector(loci=x, fitness=fitness) for x in range(n_loci)], mode=sim.ADDITIVE, begin=4, step=1)
        pre_ops = [ migrator,
                    selector,
                    sim.InitGenotype(prop=(0.1, 0.9), subPops=[2])
                    ]
        post_ops = [ sim.ResizeSubPops(subPops=[2], sizes=[math.ceil(psize*afr_size)], propagate=True, at=g) for g in range(1,11) ]

        pop = sim.Population(size=[psize,psize,math.ceil(afr_size*psize)], ploidy=2, loci=n_loci, subPopNames=['AZ', 'TX', 'AFR'], ancGen=-1, infoFields=['fitness','migrate_to', 'migrate_from']) # store all past generations

        pop.evolve(
            initOps=[sim.InitSex(maleFreq=0.9),
                     sim.InitGenotype(prop=(gene_prop_az, 1 - gene_prop_az), subPops=[0]),
                     sim.InitGenotype(prop=(gene_prop_tx, 1 - gene_prop_tx), subPops=[1]),
                     sim.InitGenotype(prop=(0.1, 0.9), subPops=[2])
                     ],
            matingScheme=matingScheme,
            preOps=[],
            postOps=[],
            gen=1
        )
        sim.stat(pop, alleleFreq=list(range(n_loci)), subPops=[0])
        az1 = np.mean(np.array([pop.dvars().alleleFreq[loc][1] for loc in range(n_loci)]))
        sim.stat(pop, alleleFreq=list(range(n_loci)), subPops=[1])
        tx1 = np.mean(np.array([pop.dvars().alleleFreq[loc][1] for loc in range(n_loci)]))

        cp = export(pop, fl, False)
        count_pop["AZ_early"] = cp["AZ"]
        count_pop["TX_early"] = cp["TX"]

        pop.evolve(
            initOps=[],
            matingScheme=matingScheme,
            preOps=pre_ops,
            postOps=post_ops,
            gen=10
        )

        sim.stat(pop, alleleFreq=list(range(n_loci)), subPops=[0])
        az2 = np.mean(np.array([pop.dvars().alleleFreq[loc][1] for loc in range(n_loci)]))
        sim.stat(pop, alleleFreq=list(range(n_loci)), subPops=[1])
        tx2 = np.mean(np.array([pop.dvars().alleleFreq[loc][1] for loc in range(n_loci)]))
        cp = export(pop, fl, True)
        count_pop["AZ_late"] = cp["AZ"]
        count_pop["TX_late"] = cp["TX"]
        print("%i)  nloci: %i      AZ : %.3f->%.3f     TX : %.3f->%.3f" % (i, n_loci, az1, az2, tx1, tx2))


    n = count_pop["AZ_early"] + count_pop["TX_early"] + count_pop["AZ_late"] + count_pop["TX_late"]
    with open(fl, 'r+') as f:
        lns = f.readlines()
        lns.insert(0, '30164 48394 29292\n')
        lns.insert(0, 'simuPOP_export %d %d\n' % (n, len(sg)))
        f.seek(0)  # readlines consumes the iterator, so we need to start over
        f.writelines(lns)

if __name__ == '__main__':
    global psize
    global afr_size
    (w, m1, m2) = (1, 0, 0)
    if args.selection:
        w = random.uniform(0.9, 1.0)
    m1 = random.uniform(0.0, 0.05)
    m2 = random.uniform(0.0, 0.2)
    psize = int(random.uniform(5000, 10000))
    afr_size = 2.0
    fl = args.ms_file
    corfl = args.corfile
    simu(w, m1, m2, psize, afr_size, fl)
    os.system("Rscript fst_sim.R %s %f %f %f %i %s" % (fl, w, m1, m2, psize, corfl))
