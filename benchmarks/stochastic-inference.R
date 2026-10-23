#!/usr/bin/env Rscript

# End-to-end stochastic generalized-EM benchmark.  The fresh bound is computed
# in a separate full E-step at the final globals for both methods, so it does
# not rely on either method's convergence-monitoring history.

suppressPackageStartupMessages(library(stm))

args <- commandArgs(trailingOnly=TRUE)
repetitions <- if(length(args)) as.integer(args[[1]]) else 1L
output <- if(length(args) >= 2) args[[2]] else NULL
if(is.na(repetitions) || repetitions < 1) stop("repetitions must be positive")

data(poliblog5k, package="stm")

prepare_case <- function(document_count, topics, prevalence=FALSE, seed=91L) {
  if(document_count == length(poliblog5k.docs)) {
    prepared <- list(documents=poliblog5k.docs, vocab=poliblog5k.voc,
                     meta=poliblog5k.meta)
  } else {
    prepared <- prepDocuments(
      poliblog5k.docs[seq_len(document_count)], poliblog5k.voc,
      poliblog5k.meta[seq_len(document_count),,drop=FALSE],
      lower.thresh=0, verbose=FALSE
    )
  }
  c(prepared, list(K=topics, prevalence=prevalence, seed=seed))
}

cases <- list(
  poliblog_1000_ctm=prepare_case(1000, 20, FALSE, 91L),
  poliblog_1000_prevalence=prepare_case(1000, 15, TRUE, 92L),
  poliblog_5000_ctm=prepare_case(5000, 30, FALSE, 91L)
)

fresh_bound <- function(model, case) {
  if(is.null(model$eta)) {
    model <- updateTheta(model, case$documents, verbose=FALSE)
  }
  timing <- system.time({
    ss <- stm:::estep(
      case$documents, model$settings$covariates$betaindex,
      !is.null(model$mu$gamma), lapply(model$beta$logbeta, exp),
      model$eta, model$mu$mu, model$sigma, FALSE
    )
  })[["elapsed"]]
  list(bound=sum(ss$bound), time=unname(timing))
}

fit_case <- function(case, method) {
  common <- list(
    documents=case$documents, vocab=case$vocab, K=case$K,
    init.type="Spectral", seed=case$seed, emtol=1e-5, verbose=FALSE
  )
  if(case$prevalence) {
    common$prevalence <- ~rating + blog
    common$data <- case$meta
  }
  arguments <- if(method == "batch") {
    c(common, list(max.em.its=120))
  } else {
    c(common, list(inference="Incremental"))
  }
  elapsed <- system.time(model <- do.call(stm, arguments))[["elapsed"]]
  score <- fresh_bound(model, case)
  data.frame(
    method=method, elapsed=unname(elapsed), fresh_bound=score$bound,
    bound_evaluation_time=score$time, iterations=model$convergence$its,
    converged=model$convergence$converged,
    updates=if(is.null(model$convergence$updates)) NA_integer_ else
      model$convergence$updates,
    batch_size=if(is.null(model$convergence$batch.size)) NA_integer_ else
      model$convergence$batch.size,
    blocks=if(is.null(model$convergence$blocks)) NA_integer_ else
      model$convergence$blocks,
    cache_mb=if(is.null(model$convergence$cache.mb)) NA_real_ else
      model$convergence$cache.mb
  )
}

results <- list()
for(case_name in names(cases)) {
  for(run in seq_len(repetitions)) {
    for(method in c("batch", "stochastic")) {
      result <- fit_case(cases[[case_name]], method)
      result$case <- case_name
      result$run <- run
      results[[length(results) + 1L]] <- result
      print(result)
    }
  }
}
results <- do.call(rbind, results)

summary <- do.call(rbind, lapply(split(results, results$case), function(x) {
  batch <- x[x$method == "batch",]
  stochastic <- x[x$method == "stochastic",]
  data.frame(
    case=x$case[1], batch_seconds=median(batch$elapsed),
    stochastic_seconds=median(stochastic$elapsed),
    speedup=median(batch$elapsed)/median(stochastic$elapsed),
    batch_bound=median(batch$fresh_bound),
    stochastic_bound=median(stochastic$fresh_bound),
    bound_delta=median(stochastic$fresh_bound)-median(batch$fresh_bound),
    stochastic_epochs=median(stochastic$iterations),
    stochastic_batch_size=median(stochastic$batch_size),
    stochastic_cache_mb=median(stochastic$cache_mb)
  )
}))
print(summary, row.names=FALSE)

if(!is.null(output)) saveRDS(list(raw=results, summary=summary), output)

