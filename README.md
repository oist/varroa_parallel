# varroa_parallel
Parallel genomic evolution of parasite tolerance by wild honey bee populations - source code.

The processing pipeline includes 3 steps: variant calling, phasing, and association.

1. Variant calling

 - Three methods are ran in parallel: freebayes, samtools, and GATK.
 - Results are next compared and filtered using BAYSIC.
 - Worker and drone data are split, variants where drones appear diploid are removed.
 - Worker data are filtered based on different amounts of missing variants.
 - Filtered worker data is merged with drone data. Before merging, synthetic 'diplod' drones are created by pairing drone individuals.

2. Phasing

 - The chomosomes are ordered by linkage group, with the correct orientation
	python split_chrom.py > commands.txt
	qsub -t 1-338 commands.sh
 - Ordered chromosomes are pooled using **merge.sh**. This script uses **running_sum.py** to re-number the positions continuously along an entire linkage group. It also uses **fix_vcf.py** to fix coding of missing genotypes, so that they are compatible with BEAGLE.
  **get_coords.py** makes a map between the official coordinates and the new coordinates
  cat coords/Group*txt > coords/all.txt
  - Haploid drones are paired one more time for testing. This also re-labels them to the new coordinate system.
	python merge_drones.py | bgzip > drones/fake_diploid.vcf.gz; tabix -p vcf drones/fake_diploid.vcf.gz
  - Phasing and imputation are performed **beagle.sh**

3. Association

  - Individuals' haplotypes are extracted for the Fst analysis **gene_haplos.py**.

