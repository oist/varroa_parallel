#!/usr/bin/env Rscript
library(PopGenome)

# Read in results of simuPOP simulations and calculate Fst values.

ns <- list("TX_early"= 98, "AZ_early" = 43, "AZ_late" = 87, "TX_late" = 102)

pop.az.early <- paste(1:(ns[["AZ_early"]]))
n <- ns[["AZ_early"]]
pop.tx.early <- paste((n+1):(n + ns[["TX_early"]]))
n <- n + ns[["TX_early"]]
pop.az.late <- paste((n+1):(n + ns[["AZ_late"]]))
n <- n + ns[["AZ_late"]]
pop.tx.late <- paste((n+1):(n + ns[["TX_late"]]))

pops <- list("az_late"=pop.az.late, "tx_late"=pop.tx.late, "az_early"=pop.az.early, "tx_early"=pop.tx.early)

comp_az <- paste("pop", which(names(pops) == "az_late"), "/pop", which(names(pops) == "az_early"), sep="")
comp_tx <- paste("pop", which(names(pops) == "tx_late"), "/pop", which(names(pops) == "tx_early"), sep="")


get_fst <- function(f, m1, m2, fl) {
  gb <- readMS(fl)

  gb <- set.populations(gb, pops)
  gb <- F_ST.stats(gb, mode="haplotype")
  fst <- gb@hap.F_ST.pairwise

  x <- fst[comp_az,]
  y <- fst[comp_tx,]
  s <- !(is.na(x) | is.na(y))
  x <- x[s]
  y <- y[s]
  cr.spearman <- cor.test(x, y, method="spearman")
  return(cr.spearman$estimate)
}

args <- commandArgs(trailingOnly=TRUE)
fl <- args[1]
f <- as.double(args[2])
m1 <- as.double(args[3])
m2 <- as.double(args[4])
psize <- as.numeric(args[5])
corfl <- args[6]

cr <- get_fst(f, m1, m2, fl)
ln <- paste(format(f, nsmall=5), format(m1,nsmall=5), format(m2,nsmall=5), psize, format(crs,nsmall=5), sep=",")

write(ln, file=corfl, append=TRUE)
