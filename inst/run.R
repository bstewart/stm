#!/usr/bin/env Rscript
library(stm)
library(stringr)

# Get command line argument for operating system
args = commandArgs(trailingOnly=TRUE)

ks <- c(5, 10, 40)
s_type <- "neum"
m_type <- "ucminf"

for(k in ks) {
  if(s_type == "neum") {
    set.seed(13)
    mod <- stm(poliblog5k.docs, poliblog5k.voc, K=k, init.type = "Spectral", max.em.its = 100, control=list(method=m_type, neum_sum_cpp=TRUE))
  }
  else {
    set.seed(13)
    mod <- stm(poliblog5k.docs, poliblog5k.voc, K=k, init.type = "Spectral", max.em.its = 100, control=list(method=m_type))
  }
  file_name = str_c(k, "_", s_type, "_", m_type, "_", args[1], ".Rds")
  saveRDS(mod, file=file_name)
}
