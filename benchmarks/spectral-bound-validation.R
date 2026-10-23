#!/usr/bin/env Rscript

# Compare final variational bounds from legacy-equivalent exact spectral
# initialization and the optimized projected-anchor initialization. Both fits
# use identical document-update code, seed, tolerance, and EM cap.

suppressPackageStartupMessages(library(stm))

args <- commandArgs(trailingOnly=TRUE)
output <- if(length(args)) args[[1]] else NULL
maximum_em_iterations <- if(length(args) >= 2) as.integer(args[[2]]) else 150L
if(is.na(maximum_em_iterations) || maximum_em_iterations < 1) {
  stop("maximum_em_iterations must be positive")
}

data(poliblog5k, package="stm")

prepare_case <- function(document_count, topics) {
  if(document_count == length(poliblog5k.docs)) {
    return(list(documents=poliblog5k.docs, vocab=poliblog5k.voc, K=topics))
  }
  prepared <- prepDocuments(
    poliblog5k.docs[seq_len(document_count)], poliblog5k.voc,
    lower.thresh=0, verbose=FALSE)
  list(documents=prepared$documents, vocab=prepared$vocab, K=topics)
}

cases <- list(
  poliblog_300=prepare_case(300, 10),
  poliblog_1000=prepare_case(1000, 20),
  poliblog_5000=prepare_case(5000, 30)
)

fit_case <- function(case, projected) {
  control <- if(projected) list() else list(SpectralRP=FALSE)
  timing <- system.time({
    fit <- stm(
      case$documents, case$vocab, K=case$K,
      max.em.its=maximum_em_iterations, emtol=1e-5,
      init.type="Spectral", control=control,
      seed=777, verbose=FALSE)
  })
  data.frame(
    method=if(projected) "optimized" else "exact_reference",
    elapsed=unname(timing[["elapsed"]]),
    iterations=fit$convergence$its,
    converged=isTRUE(fit$convergence$converged),
    bound=tail(fit$convergence$bound, 1))
}

results <- list()
for(case_name in names(cases)) {
  reference <- fit_case(cases[[case_name]], projected=FALSE)
  optimized <- fit_case(cases[[case_name]], projected=TRUE)
  result <- rbind(reference, optimized)
  result$case <- case_name
  result$D <- length(cases[[case_name]]$documents)
  result$V <- length(cases[[case_name]]$vocab)
  result$K <- cases[[case_name]]$K
  result$bound_delta <- result$bound - reference$bound
  results[[case_name]] <- result
  print(result, row.names=FALSE)
}
results <- do.call(rbind, results)
if(!is.null(output)) {
  saveRDS(results, output)
  cat("Saved bound validation to", output, "\n")
}
