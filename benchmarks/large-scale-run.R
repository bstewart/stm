#!/usr/bin/env Rscript

source(file.path("benchmarks", "large-scale-common.R"))

options <- parse_benchmark_args(list(
  method="memoized",
  corpus="large-scale-anchored-50k.rds",
  output="large-scale-result.rds",
  seed=1L,
  max_epochs=20L,
  emtol=1e-5,
  svi_defaults=FALSE,
  verbose=TRUE
))
options$method <- match.arg(
  options$method, c("memoized", "svi", "batch")
)

suppressPackageStartupMessages(library(stm))
corpus <- readRDS(options$corpus)
K <- corpus$settings$topics
common <- list(
  documents=corpus$documents,
  vocab=corpus$vocab,
  K=K,
  init.type="Custom",
  seed=options$seed,
  verbose=options$verbose,
  control=list(custom.beta=corpus$initial.logbeta)
)

cat(sprintf(
  "Starting %s: stm=%s, N=%d, K=%d, V=%d, max_epochs=%d, seed=%d\n",
  options$method, as.character(utils::packageVersion("stm")),
  length(corpus$documents), K, length(corpus$vocab),
  options$max_epochs, options$seed
))
started <- proc.time()
model <- switch(
  options$method,
  batch=do.call(stm, c(common, list(
    inference="Batch",
    max.em.its=options$max_epochs,
    emtol=options$emtol
  ))),
  memoized={
    common$control <- c(common$control, list(
      incremental.final.theta=FALSE
    ))
    do.call(stm, c(common, list(
      inference="Incremental",
      max.em.its=options$max_epochs,
      emtol=options$emtol
    )))
  },
  svi={
    if(!exists("stm_svi", mode="function")) {
      stop("The loaded stm package does not provide stm_svi(); use origin/svi")
    }
    svi.arguments <- c(
      common,
      list(
        compute_final_theta=FALSE,
        improvement_tol=options$emtol
      )
    )
    if(!options$svi_defaults) {
      svi.arguments$max_epochs <- options$max_epochs
    }
    do.call(stm_svi, svi.arguments)
  }
)
fit.seconds <- unname((proc.time() - started)[3])

cat("Running one common fresh full-corpus E-step for final comparison\n")
evaluation <- fresh_model_evaluation(
  model, corpus$documents, corpus$truth
)

summary <- data.frame(
  method=options$method,
  seed=options$seed,
  stm.version=as.character(utils::packageVersion("stm")),
  documents=length(corpus$documents),
  topics=K,
  vocabulary=length(corpus$vocab),
  tokens=corpus$tokens,
  fit.seconds=fit.seconds,
  fresh.bound=evaluation$fresh.bound,
  fresh.bound.per.token=evaluation$fresh.bound.per.token,
  evaluation.seconds=evaluation$evaluation.seconds,
  mean.hellinger=evaluation$mean.hellinger,
  max.hellinger=evaluation$max.hellinger,
  mean.anchor.recall=evaluation$mean.anchor.recall,
  minimum.anchor.recall=evaluation$minimum.anchor.recall,
  mean.theta.correlation=evaluation$mean.theta.correlation,
  minimum.theta.correlation=evaluation$minimum.theta.correlation,
  iterations=if(is.null(model$convergence$its)) NA_real_ else
    model$convergence$its,
  epochs=if(!is.null(model$convergence$epochs)) {
    model$convergence$epochs
  } else if(!is.null(model$convergence$epoch)) {
    model$convergence$epoch
  } else {
    NA_real_
  },
  updates=if(is.null(model$convergence$updates)) NA_real_ else
    model$convergence$updates,
  converged=isTRUE(model$convergence$converged),
  stringsAsFactors=FALSE
)

model$theta <- NULL
model$eta <- NULL
result <- list(
  summary=summary,
  convergence=model$convergence,
  beta=model$beta,
  mu=model$mu,
  sigma=model$sigma,
  assignment=evaluation$assignment,
  benchmark.options=options,
  corpus.settings=corpus$settings,
  model.settings=model$settings
)
saveRDS(result, options$output, compress=FALSE)
utils::write.csv(summary, sub("[.]rds$", ".csv", options$output),
                 row.names=FALSE)
print(summary)

