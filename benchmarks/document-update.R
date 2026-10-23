#!/usr/bin/env Rscript

# Reproducible microbenchmark for STM's document-level variational update.
#
# Run from the package root with an R library containing RcppArmadillo:
#   R_LIBS=/path/to/library Rscript benchmarks/document-update.R

suppressPackageStartupMessages({
  library(Rcpp)
  library(RcppArmadillo)
})

args <- commandArgs(trailingOnly = TRUE)
quick <- "--quick" %in% args
quick_k100 <- "--quick-k100" %in% args
profile <- "--profile" %in% args
fused_reltol <- suppressWarnings(as.numeric(Sys.getenv("STM_FUSED_RELTOL", NA)))

if (!file.exists("src/STMCfuns.cpp") || !file.exists("data/poliblog5k.RData")) {
  stop("Run this benchmark from the stm package root.")
}

Rcpp::sourceCpp("src/STMCfuns.cpp", rebuild = TRUE, showOutput = FALSE)
load("data/poliblog5k.RData")

make_problem <- function(K, n_docs, seed) {
  set.seed(seed)
  V <- length(poliblog5k.voc)
  document_ids <- sample.int(length(poliblog5k.docs), n_docs)
  documents <- poliblog5k.docs[document_ids]

  beta <- matrix(rgamma(K * V, shape = 0.1), nrow = K, ncol = V)
  beta <- beta / rowSums(beta)

  rho <- 0.35
  sigma <- rho ^ abs(outer(seq_len(K - 1), seq_len(K - 1), `-`))
  sigma <- sigma + diag(0.25, K - 1)
  sigma_chol <- chol(sigma)

  list(
    K = K,
    documents = documents,
    beta = beta,
    mu = rnorm(K - 1, sd = 0.15),
    siginv = chol2inv(sigma_chol),
    sigmaentropy = sum(log(diag(sigma_chol))),
    eta = rep(0, K - 1),
    start = "cold"
  )
}

available_engines <- function() {
  engines <- "legacy"
  if (exists("gradcpp_nocopy", mode = "function")) {
    engines <- c(engines, "nocopy")
  }
  if (exists("hpbcpp_potri", mode = "function")) {
    engines <- c(engines, "potri")
  }
  if (exists("optimizecpp", mode = "function")) {
    engines <- c(engines, "fused")
  }
  engines
}

optimize_one <- function(problem, document, eta, engine, control) {
  words <- document[1, ]
  counts <- document[2, ]
  beta_document <- problem$beta[, words, drop = FALSE]

  if (engine == "fused") {
    if (!is.na(fused_reltol)) control$reltol <- fused_reltol
    return(optimizecpp(
      eta = eta,
      beta = beta_document,
      doc_ct = counts,
      mu = problem$mu,
      siginv = problem$siginv,
      control = control
    ))
  }

  gradient <- if (engine %in% c("nocopy", "potri")) {
    gradcpp_nocopy
  } else {
    gradcpp
  }
  fit <- optim(
    par = eta,
    fn = lhoodcpp,
    gr = gradient,
    method = "BFGS",
    control = control,
    beta = beta_document,
    doc_ct = counts,
    mu = problem$mu,
    siginv = problem$siginv
  )
  list(
    par = fit$par,
    value = fit$value,
    counts = unname(fit$counts),
    convergence = fit$convergence
  )
}

posterior_one <- function(problem, document, eta, engine) {
  words <- document[1, ]
  beta_document <- problem$beta[, words, drop = FALSE]
  posterior <- if (engine %in% c("potri", "fused") &&
                   exists("hpbcpp_potri", mode = "function")) {
    hpbcpp_potri
  } else if (engine != "legacy" &&
             exists("hpbcpp_nocopy", mode = "function")) {
    hpbcpp_nocopy
  } else {
    hpbcpp
  }
  posterior(
    eta = eta,
    beta = beta_document,
    doc_ct = document[2, ],
    mu = problem$mu,
    siginv = problem$siginv,
    sigmaentropy = problem$sigmaentropy
  )
}

run_engine <- function(problem, engine, control) {
  fits <- vector("list", length(problem$documents))
  optimization_time <- system.time({
    for (i in seq_along(problem$documents)) {
      initial_eta <- if (is.list(problem$eta)) problem$eta[[i]] else problem$eta
      fits[[i]] <- optimize_one(
        problem, problem$documents[[i]], initial_eta, engine, control
      )
    }
  })[["elapsed"]]

  posteriors <- vector("list", length(problem$documents))
  posterior_time <- system.time({
    for (i in seq_along(problem$documents)) {
      posteriors[[i]] <- posterior_one(
        problem, problem$documents[[i]], fits[[i]]$par, engine
      )
    }
  })[["elapsed"]]

  objectives <- numeric(length(fits))
  gradient_norms <- numeric(length(fits))
  for (i in seq_along(fits)) {
    document <- problem$documents[[i]]
    beta_document <- problem$beta[, document[1, ], drop = FALSE]
    objectives[i] <- lhoodcpp(
      fits[[i]]$par, beta_document, document[2, ],
      problem$mu, problem$siginv
    )
    gradient_norms[i] <- max(abs(gradcpp(
      fits[[i]]$par, beta_document, document[2, ],
      problem$mu, problem$siginv
    )))
  }

  counts <- do.call(rbind, lapply(fits, `[[`, "counts"))
  data.frame(
    engine = engine,
    start = problem$start,
    K = problem$K,
    documents = length(problem$documents),
    median_unique_terms = median(vapply(problem$documents, ncol, integer(1))),
    optimization_seconds = optimization_time,
    posterior_seconds = posterior_time,
    total_seconds = optimization_time + posterior_time,
    objective_sum = sum(objectives),
    bound_sum = sum(vapply(posteriors, `[[`, numeric(1), "bound")),
    max_gradient_inf_norm = max(gradient_norms),
    function_evaluations = sum(counts[, 1]),
    gradient_evaluations = sum(counts[, 2]),
    failures = sum(vapply(fits, `[[`, numeric(1), "convergence") != 0),
    stringsAsFactors = FALSE
  )
}

make_warm_problem <- function(problem, control, seed) {
  message("Preparing warm starts for K=", problem$K)
  modes <- vector("list", length(problem$documents))
  for (i in seq_along(problem$documents)) {
    modes[[i]] <- optimize_one(
      problem, problem$documents[[i]], problem$eta, "legacy", control
    )$par
  }
  set.seed(seed)
  problem$eta <- lapply(modes, function(mode) {
    mode + rnorm(length(mode), sd = 0.02)
  })
  problem$start <- "warm"
  problem
}

profile_legacy <- function(problem, control) {
  output <- tempfile("stm-document-update-profile-", fileext = ".out")
  on.exit(unlink(output), add = TRUE)
  Rprof(output, interval = 0.0005)
  on.exit(Rprof(NULL), add = TRUE)
  invisible(run_engine(problem, "legacy", control))
  Rprof(NULL)
  summaryRprof(output)
}

control <- list(maxit = 500)
cases <- if (quick_k100) {
  list(make_problem(K = 100, n_docs = 200, seed = 7103))
} else if (quick) {
  list(make_problem(K = 20, n_docs = 40, seed = 7101))
} else {
  cold_cases <- list(
    make_problem(K = 20, n_docs = 5000, seed = 7101),
    make_problem(K = 50, n_docs = 3000, seed = 7102),
    make_problem(K = 100, n_docs = 1500, seed = 7103)
  )
  warm_cases <- Map(
    function(problem, seed) make_warm_problem(problem, control, seed),
    cold_cases, c(8101, 8102, 8103)
  )
  c(cold_cases, warm_cases)
}
engines <- available_engines()

# Warm compiled code and BLAS before timing.
invisible(optimize_one(cases[[1]], cases[[1]]$documents[[1]], cases[[1]]$eta,
                       engines[[1]], control))

n_repetitions <- if (quick || quick_k100) 1L else 3L
results <- do.call(rbind, lapply(seq_len(n_repetitions), function(repetition) {
  engine_order <- if (repetition %% 2L) engines else rev(engines)
  repetition_results <- do.call(rbind, lapply(cases, function(problem) {
    do.call(rbind, lapply(engine_order, function(engine) {
      message("repetition=", repetition, ", K=", problem$K,
              ", engine=", engine)
      result <- run_engine(problem, engine, control)
      result$repetition <- repetition
      result
    }))
  }))
  repetition_results
}))

group_key <- interaction(results$K, results$start, results$repetition,
                         drop = TRUE)
legacy_rows <- results$engine == "legacy"
legacy_objective <- setNames(results$objective_sum[legacy_rows],
                             group_key[legacy_rows])
legacy_bound <- setNames(results$bound_sum[legacy_rows], group_key[legacy_rows])
results$objective_delta <- results$objective_sum - legacy_objective[group_key]
results$bound_delta <- results$bound_sum - legacy_bound[group_key]

print(results, row.names = FALSE, digits = 12)

if (!quick && !quick_k100) {
  medians <- aggregate(
    cbind(optimization_seconds, posterior_seconds, total_seconds) ~
      engine + K + start,
    data = results,
    FUN = median
  )
  baseline <- setNames(
    medians$total_seconds[medians$engine == "legacy"],
    interaction(medians$K[medians$engine == "legacy"],
                medians$start[medians$engine == "legacy"])
  )
  medians$speedup_vs_legacy <-
    unname(baseline[interaction(medians$K, medians$start)]) /
    medians$total_seconds
  message("Median timings")
  print(medians[order(medians$K, medians$start, medians$engine), ],
        row.names = FALSE, digits = 6)
}

if (profile) {
  prof <- profile_legacy(cases[[min(2, length(cases))]], control)
  print(head(prof$by.total, 20))
  print(head(prof$by.self, 20))
}
