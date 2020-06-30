#!/usr/bin/env Rscript
library(stm)
library(stringr)

# Get command line argument for operating system
args = commandArgs(trailingOnly=TRUE)

# Generate data
docs = poliblog5k.docs
voc = poliblog5k.voc
heldout <- make.heldout(docs, voc)

reg <- stm(heldout$documents, heldout$vocab, K=5, init.type = "Spectral", max.em.its = 100, control=list(method="BFGS"))

reg_file = str_c("reg_", args[1], ".Rda")
save(reg, file=reg_file)

neum <- stm(heldout$documents, heldout$vocab, K=5, init.type = "Spectral", max.em.its = 100, control=list(method="BFGS", neum_sum_cpp=TRUE))

neum_file = str_c("neum_", args[1], ".Rda")
save(neum, file=neum_file)

