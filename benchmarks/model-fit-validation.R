#!/usr/bin/env Rscript

# End-to-end numerical validation for the document update. Run this script in
# separate R processes with the baseline and candidate package libraries.

suppressPackageStartupMessages(library(stm))
args <- commandArgs(trailingOnly = TRUE)
output <- if (length(args)) args[[1]] else tempfile(fileext = ".rds")

data(poliblog5k, package = "stm")
documents <- poliblog5k.docs[seq_len(1000)]

timing <- system.time({
  fit <- stm(
    documents = documents,
    vocab = poliblog5k.voc,
    K = 20,
    max.em.its = 12,
    init.type = "Random",
    seed = 9017,
    verbose = FALSE
  )
})

result <- list(
  elapsed = unname(timing[["elapsed"]]),
  bound = fit$convergence$bound,
  eta = fit$eta,
  sigma = fit$sigma,
  beta = fit$beta$logbeta
)
saveRDS(result, output)
print(list(elapsed = result$elapsed, bound = result$bound))
