#!/usr/bin/env Rscript
library(stm)

# Get command line arguments 1 doc order, 2 summation
args = commandArgs(trailingOnly=TRUE)

# Load data set settings
prepped.docs <- readRDS("prepped_docs.rds")
num_topics=105
max_its=25

set.seed(08540)

if(args[1]=="rand" & args[2]=="neum") {
  cat("Running rand/neum... \n")
  rand_neum <- stm(prepped.docs$documents,
                   prepped.docs$vocab,
                   K=num_topics,
                   max.em.its = max_its,
                   prevalence = ~ matched + s(date),
                   data = prepped.docs$meta,
                   init.type = "Spectral",
                   verbose = FALSE,
                   control=list(method="BFGS", rand_docs=TRUE, neum_sum_cpp=TRUE))
  save(rand_neum,file="rand_neum.Rda")
}

if(args[1]=="deter" & args[2]=="neum") {
  cat("Running deter/neum... \n")
  deter_neum <- stm(prepped.docs$documents,
                   prepped.docs$vocab,
                   K=num_topics,
                   max.em.its = max_its,
                   prevalence = ~ matched + s(date),
                   data = prepped.docs$meta,
                   init.type = "Spectral",
                   verbose = FALSE,
                   control=list(method="BFGS", neum_sum_cpp=TRUE))
  save(deter_neum,file="deter_neum.Rda")
}

if(args[1]=="rand" & args[2]=="reg") {
  cat("Running rand/reg... \n")
  rand_reg <- stm(prepped.docs$documents,
                    prepped.docs$vocab,
                    K=num_topics,
                    max.em.its = max_its,
                    prevalence = ~ matched + s(date),
                    data = prepped.docs$meta,
                    init.type = "Spectral",
                    verbose = FALSE,
                    control=list(method="BFGS", rand_docs=TRUE))
  save(rand_reg,file="rand_reg.Rda")
}

if(args[1]=="deter" & args[2]=="reg") {
  cat("Running deter/reg... \n")
  deter_reg <- stm(prepped.docs$documents,
                    prepped.docs$vocab,
                    K=num_topics,
                    max.em.its = max_its,
                    prevalence = ~ matched + s(date),
                    data = prepped.docs$meta,
                    init.type = "Spectral",
                    verbose = FALSE,
                    control=list(method="BFGS"))
  save(deter_reg,file="deter_reg.Rda")
}