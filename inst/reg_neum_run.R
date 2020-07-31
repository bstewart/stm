#!/usr/bin/env Rscript
library(stm)
library(stringr)

# Get command line argument for operating system
args = commandArgs(trailingOnly=TRUE)

set.seed(13)
reg <- stm(poliblog5k.docs, poliblog5k.voc, K=5, init.type = "Spectral", max.em.its = 100, control=list(method="ucminf"))
reg_file = str_c("reg_", args[1], ".Rds")
saveRDS(reg, file=reg_file)

set.seed(13)
neum <- stm(poliblog5k.docs, poliblog5k.voc, K=5, init.type = "Spectral", max.em.its = 100, control=list(method="ucminf", neum_sum_cpp=TRUE))
neum_file = str_c("neum_", args[1], ".Rds")
saveRDS(neum, file=neum_file)

