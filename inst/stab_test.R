#!/usr/bin/env Rscript
library(stm)
library(stringr)

# Get command line argument for operating system
args = commandArgs(trailingOnly=TRUE)

reg <- stm(poliblog5k.docs, poliblog5k.voc, K=5, init.type = "Spectral", max.em.its = 100, control=list(method="BFGS"))

reg_file = str_c("reg_", args[1], ".Rda")
save(reg, file=reg_file)

neum <- stm(poliblog5k.docs, poliblog5k.voc, K=5, init.type = "Spectral", max.em.its = 100, control=list(method="BFGS", neum_sum_cpp=TRUE))

neum_file = str_c("neum_", args[1], ".Rda")
save(neum, file=neum_file)

