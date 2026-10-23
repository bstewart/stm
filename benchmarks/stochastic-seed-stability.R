#!/usr/bin/env Rscript

# Assess stochastic-inference stability separately for
#   1. the stochastic block/order seed,
#   2. the randomized spectral-projection seed, and
#   3. paired changes to both seeds.
#
# The projected spectral initializer is randomized in construction but
# deterministic by default: rp.seed defaults to 9 and is independent of the
# model seed.  This benchmark changes rp.seed explicitly when sketch variation
# is desired.
#
# Usage:
#   Rscript benchmarks/stochastic-seed-stability.R [nseeds] [output.rds]

suppressPackageStartupMessages(library(stm))

args <- commandArgs(trailingOnly = TRUE)
nseeds <- if (length(args)) as.integer(args[[1]]) else 5L
if (is.na(nseeds) || nseeds < 3L || nseeds > 5L) {
  stop("nseeds must be between 3 and 5")
}
output <- if (length(args) > 1L) args[[2]] else
  file.path(tempdir(), "stm-stochastic-seed-stability.rds")

model_seed <- 91L
default_rp_seed <- 9L
default_stochastic_seed <- model_seed + 104729L
rp_seeds <- c(default_rp_seed, 17L, 41L, 89L, 151L)[seq_len(nseeds)]
stochastic_seeds <- c(default_stochastic_seed, 2003L, 3001L, 4001L,
                      5003L)[seq_len(nseeds)]

fresh_bound <- function(model, documents) {
  if (is.null(model$eta)) {
    model <- updateTheta(model, documents, verbose = FALSE)
  }
  fresh <- stm:::estep(
    documents,
    model$settings$covariates$betaindex,
    !is.null(model$mu$gamma),
    lapply(model$beta$logbeta, exp),
    model$eta,
    model$mu$mu,
    model$sigma,
    FALSE
  )
  sum(fresh$bound)
}

# Shortest augmenting-path assignment solver.  This keeps the benchmark
# self-contained when the Suggested package clue is unavailable.
solve_assignment <- function(cost) {
  cost <- as.matrix(cost)
  n <- nrow(cost)
  m <- ncol(cost)
  stopifnot(n > 0L, n <= m, all(is.finite(cost)))
  u <- numeric(n + 1L)
  v <- numeric(m + 1L)
  p <- integer(m + 1L)
  way <- integer(m + 1L)

  for (i in seq_len(n)) {
    p[[1L]] <- i
    j0 <- 1L
    minv <- rep(Inf, m + 1L)
    used <- rep(FALSE, m + 1L)
    repeat {
      used[[j0]] <- TRUE
      i0 <- p[[j0]]
      delta <- Inf
      j1 <- 1L
      for (j in 2L:(m + 1L)) {
        if (!used[[j]]) {
          cur <- cost[i0, j - 1L] - u[[i0 + 1L]] - v[[j]]
          if (cur < minv[[j]]) {
            minv[[j]] <- cur
            way[[j]] <- j0
          }
          if (minv[[j]] < delta) {
            delta <- minv[[j]]
            j1 <- j
          }
        }
      }
      for (j in seq_len(m + 1L)) {
        if (used[[j]]) {
          u[[p[[j]] + 1L]] <- u[[p[[j]] + 1L]] + delta
          v[[j]] <- v[[j]] - delta
        } else {
          minv[[j]] <- minv[[j]] - delta
        }
      }
      j0 <- j1
      if (p[[j0]] == 0L) break
    }
    repeat {
      j1 <- way[[j0]]
      p[[j0]] <- p[[j1]]
      j0 <- j1
      if (j0 == 1L) break
    }
  }

  assignment <- integer(n)
  for (j in 2L:(m + 1L)) {
    if (p[[j]] > 0L) assignment[[p[[j]]]] <- j - 1L
  }
  assignment
}

# Validate the local assignment solver against exhaustive search on a small
# deterministic problem before it is used for topic matching.
all_permutations <- function(x) {
  if (length(x) == 1L) return(matrix(x, nrow = 1L))
  do.call(rbind, lapply(seq_along(x), function(i) {
    cbind(x[[i]], all_permutations(x[-i]))
  }))
}
set.seed(101)
assignment_check <- matrix(runif(25L), 5L, 5L)
assignment <- solve_assignment(-assignment_check)
permutations <- all_permutations(seq_len(5L))
exhaustive <- apply(permutations, 1L, function(x) {
  sum(assignment_check[cbind(seq_len(5L), x)])
})
stopifnot(abs(sum(assignment_check[cbind(seq_len(5L), assignment)]) -
              max(exhaustive)) < 1e-12)

topic_metrics <- function(reference, candidate, top_n = 20L) {
  reference <- as.matrix(reference)
  candidate <- as.matrix(candidate)
  affinity <- sqrt(reference) %*% t(sqrt(candidate))
  assignment <- solve_assignment(-affinity)
  matched <- candidate[assignment, , drop = FALSE]
  matched_affinity <- rowSums(sqrt(reference * matched))
  hellinger <- sqrt(pmax(0, 1 - matched_affinity))
  total_variation <- 0.5 * rowSums(abs(reference - matched))
  top_overlap <- vapply(seq_len(nrow(reference)), function(k) {
    ref_top <- order(reference[k, ], decreasing = TRUE)[seq_len(top_n)]
    candidate_top <- order(matched[k, ], decreasing = TRUE)[seq_len(top_n)]
    length(intersect(ref_top, candidate_top))
  }, numeric(1L))
  list(
    assignment = assignment,
    mean_affinity = mean(matched_affinity),
    min_affinity = min(matched_affinity),
    mean_hellinger = mean(hellinger),
    max_hellinger = max(hellinger),
    mean_total_variation = mean(total_variation),
    max_total_variation = max(total_variation),
    mean_top20_overlap = mean(top_overlap),
    min_top20_overlap = min(top_overlap)
  )
}

theta_metrics <- function(reference, candidate, assignment) {
  correlations <- vapply(seq_along(assignment), function(k) {
    suppressWarnings(cor(reference[, k], candidate[, assignment[[k]]]))
  }, numeric(1L))
  correlations <- correlations[is.finite(correlations)]
  c(mean_theta_correlation = mean(correlations),
    min_theta_correlation = min(correlations))
}

make_case <- function(name, n, k, prevalence = FALSE) {
  data(poliblog5k, package = "stm")
  documents <- poliblog5k.docs[seq_len(n)]
  meta <- poliblog5k.meta[seq_len(n), , drop = FALSE]
  list(
    name = name,
    documents = documents,
    vocab = poliblog5k.voc,
    data = meta,
    K = k,
    prevalence = if (prevalence) ~ rating + blog else NULL,
    ntokens = sum(vapply(documents, function(x) sum(x[2L, ]), numeric(1L)))
  )
}

cases <- list(
  make_case("D1000-K20-CTM", 1000L, 20L),
  make_case("D1000-K15-prevalence", 1000L, 15L, prevalence = TRUE),
  make_case("D5000-K30-CTM", 5000L, 30L)
)
names(cases) <- vapply(cases, `[[`, character(1L), "name")

fit_cache <- new.env(parent = emptyenv())
fit_case <- function(case, rp_seed, stochastic_seed) {
  key <- paste(case$name, rp_seed, stochastic_seed, sep = "-")
  if (exists(key, fit_cache, inherits = FALSE)) {
    return(get(key, fit_cache, inherits = FALSE))
  }
  message(sprintf("Fitting %-22s rp.seed=%-4d incremental.seed=%d",
                  case$name, rp_seed, stochastic_seed))
  fit_args <- list(
    documents = case$documents,
    vocab = case$vocab,
    K = case$K,
    data = case$data,
    init.type = "Spectral",
    inference = "Incremental",
    seed = model_seed,
    max.em.its = 75L,
    control = list(
      rp.seed = as.integer(rp_seed),
      incremental.seed = as.integer(stochastic_seed),
      incremental.final.theta = TRUE
    ),
    verbose = FALSE
  )
  if (!is.null(case$prevalence)) fit_args$prevalence <- case$prevalence
  elapsed <- system.time(model <- do.call(stm, fit_args))[["elapsed"]]
  result <- list(
    model = model,
    beta = exp(model$beta$logbeta[[1L]]),
    theta = model$theta,
    bound = fresh_bound(model, case$documents),
    elapsed = unname(elapsed),
    epochs = model$convergence$its,
    converged = isTRUE(model$convergence$converged),
    rp_seed = rp_seed,
    stochastic_seed = stochastic_seed
  )
  assign(key, result, fit_cache)
  result
}

compare_group <- function(fits, case, source) {
  reference <- fits[[1L]]
  do.call(rbind, lapply(seq_along(fits), function(i) {
    fit <- fits[[i]]
    topics <- topic_metrics(reference$beta, fit$beta)
    theta <- theta_metrics(reference$theta, fit$theta, topics$assignment)
    data.frame(
      case = case$name,
      source = source,
      seed_index = i,
      rp_seed = fit$rp_seed,
      stochastic_seed = fit$stochastic_seed,
      elapsed = fit$elapsed,
      epochs = fit$epochs,
      converged = fit$converged,
      bound = fit$bound,
      bound_per_token = fit$bound / case$ntokens,
      mean_affinity = topics$mean_affinity,
      min_affinity = topics$min_affinity,
      mean_hellinger = topics$mean_hellinger,
      max_hellinger = topics$max_hellinger,
      mean_total_variation = topics$mean_total_variation,
      max_total_variation = topics$max_total_variation,
      mean_top20_overlap = topics$mean_top20_overlap,
      min_top20_overlap = topics$min_top20_overlap,
      mean_theta_correlation = theta[["mean_theta_correlation"]],
      min_theta_correlation = theta[["min_theta_correlation"]]
    )
  }))
}

summarize_group <- function(x) {
  comparisons <- x[x$seed_index > 1L, , drop = FALSE]
  data.frame(
    case = x$case[[1L]],
    source = x$source[[1L]],
    nseeds = nrow(x),
    converged = sum(x$converged),
    bound_mean = mean(x$bound),
    bound_sd = sd(x$bound),
    bound_range = diff(range(x$bound)),
    bound_per_token_range = diff(range(x$bound_per_token)),
    elapsed_mean = mean(x$elapsed),
    epoch_min = min(x$epochs),
    epoch_max = max(x$epochs),
    mean_hellinger = mean(comparisons$mean_hellinger),
    worst_topic_hellinger = max(comparisons$max_hellinger),
    mean_total_variation = mean(comparisons$mean_total_variation),
    worst_topic_total_variation = max(comparisons$max_total_variation),
    mean_top20_overlap = mean(comparisons$mean_top20_overlap),
    worst_topic_top20_overlap = min(comparisons$min_top20_overlap),
    mean_theta_correlation = mean(comparisons$mean_theta_correlation),
    worst_topic_theta_correlation = min(comparisons$min_theta_correlation)
  )
}

# Directly test the point that the default projected initializer does not vary
# with the ordinary model seed.  It changes only when rp.seed changes.
initial_beta <- function(case, seed, rp_seed) {
  fit_args <- list(
    documents = case$documents,
    vocab = case$vocab,
    K = case$K,
    data = case$data,
    init.type = "Spectral",
    inference = "Batch",
    seed = as.integer(seed),
    max.em.its = 0L,
    control = list(rp.seed = as.integer(rp_seed)),
    verbose = FALSE
  )
  if (!is.null(case$prevalence)) fit_args$prevalence <- case$prevalence
  exp(do.call(stm, fit_args)$beta$logbeta[[1L]])
}

small_case <- cases[["D1000-K20-CTM"]]
fixed_rp_initial <- lapply(c(1L, model_seed, 999L), function(seed) {
  initial_beta(small_case, seed, default_rp_seed)
})
default_initializer_max_abs_difference <- max(vapply(
  fixed_rp_initial[-1L],
  function(x) max(abs(x - fixed_rp_initial[[1L]])),
  numeric(1L)
))
varying_rp_initial <- lapply(rp_seeds, function(seed) {
  initial_beta(small_case, model_seed, seed)
})
initialization_stability <- do.call(rbind, lapply(seq_along(rp_seeds), function(i) {
  metrics <- topic_metrics(varying_rp_initial[[1L]], varying_rp_initial[[i]])
  data.frame(
    seed_index = i,
    rp_seed = rp_seeds[[i]],
    mean_hellinger = metrics$mean_hellinger,
    max_hellinger = metrics$max_hellinger,
    mean_total_variation = metrics$mean_total_variation,
    mean_top20_overlap = metrics$mean_top20_overlap,
    min_top20_overlap = metrics$min_top20_overlap
  )
}))

# Decompose the two randomness sources on the medium CTM case.
schedule_fits <- lapply(stochastic_seeds, function(seed) {
  fit_case(small_case, default_rp_seed, seed)
})
schedule_results <- compare_group(schedule_fits, small_case, "update-seed-only")

sketch_fits <- lapply(rp_seeds, function(seed) {
  fit_case(small_case, seed, default_stochastic_seed)
})
sketch_results <- compare_group(sketch_fits, small_case, "sketch-seed-only")

# Change both sources together across all three benchmark cases.
paired_results <- do.call(rbind, lapply(cases, function(case) {
  fits <- lapply(seq_len(nseeds), function(i) {
    fit_case(case, rp_seeds[[i]], stochastic_seeds[[i]])
  })
  compare_group(fits, case, "paired-seeds")
}))

raw <- rbind(schedule_results, sketch_results, paired_results)
groups <- split(raw, interaction(raw$case, raw$source, drop = TRUE))
summary <- do.call(rbind, lapply(groups, summarize_group))
rownames(summary) <- NULL

result <- list(
  generated = Sys.time(),
  stm_version = as.character(packageVersion("stm")),
  model_seed = model_seed,
  rp_seeds = rp_seeds,
  stochastic_seeds = stochastic_seeds,
  default_initializer_max_abs_difference =
    default_initializer_max_abs_difference,
  initialization_stability = initialization_stability,
  raw = raw,
  summary = summary
)
saveRDS(result, output)

cat("\nDefault fixed-sketch initializer max absolute difference across model seeds:",
    format(default_initializer_max_abs_difference, scientific = TRUE), "\n\n")
cat("Initialization stability when rp.seed changes:\n")
print(initialization_stability, row.names = FALSE, digits = 4)
cat("\nFinal-fit seed stability:\n")
print(summary, row.names = FALSE, digits = 5)
cat("\nSaved raw results to", normalizePath(output, mustWork = FALSE), "\n")
