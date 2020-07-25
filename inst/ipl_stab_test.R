#!/usr/bin/env Rscript
library(stm)
library(stringr)

# Get command line argument for operating system
args = commandArgs(trailingOnly=TRUE)

cat(paste("OS:", args[1]))

set.seed(13)
reg <- stm(poliblog5k.docs, poliblog5k.voc, K=5, init.type = "Spectral", seed = 13, max.em.its = 100, control=list(method="BFGS"))
reg_rand <- stm(poliblog5k.docs, poliblog5k.voc, K=5, init.type = "Spectral", seed = 13, max.em.its = 100, control=list(method="BFGS", rand_docs=TRUE)) 
cat(paste("Compare regular summation:\n"))
print(all.equal(reg, reg_rand))

set.seed(13)
neum <- stm(poliblog5k.docs, poliblog5k.voc, K=5, init.type = "Spectral", seed = 13, max.em.its = 100, control=list(method="BFGS", neum_sum_cpp=TRUE))
neum_rand <- stm(poliblog5k.docs, poliblog5k.voc, K=5, init.type = "Spectral", seed = 13, max.em.its = 100, control=list(method="BFGS", rand_docs=TRUE, neum_sum_cpp=TRUE)) 
cat(paste("Compare Neumaier summation:\n"))
print(all.equal(neum, neum_rand))
