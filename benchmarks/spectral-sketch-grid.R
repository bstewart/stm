#!/usr/bin/env Rscript

# Width/density ablation for the randomized spectral anchor sketch.
#
# The projected representation is used only for anchor selection.  Every
# configuration recovers beta from exact original-space kernels, then runs the
# same stochastic optimizer with a fixed block seed.  This separates sketch
# geometry from downstream stochastic-update randomness.
#
# The D=1,000 case receives the full grid.  The larger D=5,000 case receives
# four prespecified corner configurations: default, width only, density only,
# and width plus density.
#
# Usage:
#   Rscript benchmarks/spectral-sketch-grid.R [nseeds] [output.rds] [mode]
#                                             [full-grid.rds]
#
# mode is "full" (default) or "large-middle" (D=5,000 with p512-s05 and
# p512-s25 only).

suppressPackageStartupMessages(library(stm))

args <- commandArgs(trailingOnly = TRUE)
nseeds <- if (length(args)) as.integer(args[[1L]]) else 8L
if (is.na(nseeds) || nseeds < 5L || nseeds > 8L) {
  stop("nseeds must be between 5 and 8")
}
output <- if (length(args) > 1L) args[[2L]] else
  file.path(tempdir(), "stm-spectral-sketch-grid.rds")
mode <- if (length(args) > 2L) args[[3L]] else "full"
if (!mode %in% c("full", "large-middle")) {
  stop("mode must be 'full' or 'large-middle'")
}
baseline_file <- if (length(args) > 3L) args[[4L]] else NULL
if (mode == "large-middle" && is.null(baseline_file)) {
  stop("large-middle mode requires the completed full-grid RDS as argument 4")
}

rp_seeds <- c(9L, 17L, 41L, 89L, 151L, 3L, 5L, 7L)[seq_len(nseeds)]
model_seed <- 91L
stochastic_seed <- model_seed + 104729L

configurations <- data.frame(
  config = c(
    "p256-s05", "p512-s05", "p1024-s05",
    "p256-s25", "p512-s25", "p1024-s25", "p256-dense"
  ),
  p = c(256L, 512L, 1024L, 256L, 512L, 1024L, 256L),
  density = c(.05, .05, .05, .25, .25, .25, 1),
  large_case = c(TRUE, TRUE, TRUE, FALSE, TRUE, TRUE, TRUE),
  stringsAsFactors = FALSE
)
configurations$nonzeros <- with(
  configurations, pmin(p, pmax(1L, ceiling(p * density)))
)

make_case <- function(name, n, k) {
  data(poliblog5k, package = "stm")
  documents <- poliblog5k.docs[seq_len(n)]
  list(
    name = name,
    documents = documents,
    vocab = poliblog5k.voc,
    K = k,
    ntokens = sum(vapply(documents, function(x) sum(x[2L, ]), numeric(1L)))
  )
}

cases <- list(
  make_case("D1000-K20-CTM", 1000L, 20L),
  make_case("D5000-K30-CTM", 5000L, 30L)
)
names(cases) <- vapply(cases, `[[`, character(1L), "name")
if (mode == "large-middle") {
  cases <- cases["D5000-K30-CTM"]
  configurations <- configurations[
    configurations$config %in% c("p512-s05", "p512-s25"), , drop = FALSE
  ]
}

restore_words <- function(beta, valid, vocabulary_size) {
  if (all(valid)) return(beta)
  if (ncol(beta) == sum(valid)) {
    restored <- matrix(0, nrow = nrow(beta), ncol = vocabulary_size)
    restored[, valid] <- beta
  } else {
    restored <- beta
  }
  restored <- restored + .001 / vocabulary_size
  restored / rowSums(restored)
}

projected_initialization <- function(case, p, density, seed) {
  elapsed <- system.time({
    sketch <- stm:::spectralProjectionCpp(
      case$documents, length(case$vocab),
      projection_dimension = as.integer(p),
      projection_sparsity = density,
      seed = as.integer(seed)
    )
    valid <- sketch$row_normalizer > 0
    anchor_reduced <- stm:::fastAnchorCpp(
      sketch$Qbar[valid, , drop = FALSE], case$K
    )
    anchors <- which(valid)[anchor_reduced]
    beta <- stm:::recoverL2OriginalCpp(
      case$documents, length(case$vocab), anchors,
      sketch$row_normalizer, sketch$wprob,
      tolerance = 1e-6, maximum_iterations = 500L
    )$A
    beta <- restore_words(beta, valid, length(case$vocab))
  })[["elapsed"]]
  list(beta = beta, anchors = anchors, elapsed = unname(elapsed))
}

exact_initialization <- function(case) {
  elapsed <- system.time({
    moment <- stm:::spectralGramCpp(case$documents, length(case$vocab))
    valid <- rowSums(moment$Q) > 0
    qbar <- moment$Q[valid, valid, drop = FALSE] / rowSums(moment$Q)[valid]
    anchor_reduced <- stm:::fastAnchorCpp(qbar, case$K)
    anchors <- which(valid)[anchor_reduced]
    beta <- stm:::recoverL2Cpp(
      qbar, anchor_reduced, moment$wprob[valid],
      tolerance = 1e-7, maximum_iterations = 500L
    )$A
    beta <- restore_words(beta, valid, length(case$vocab))
  })[["elapsed"]]
  list(beta = beta, anchors = anchors, elapsed = unname(elapsed))
}

fresh_bound <- function(model, documents) {
  if (is.null(model$eta)) {
    model <- updateTheta(model, documents, verbose = FALSE)
  }
  sum(stm:::estep(
    documents,
    model$settings$covariates$betaindex,
    !is.null(model$mu$gamma),
    lapply(model$beta$logbeta, exp),
    model$eta,
    model$mu$mu,
    model$sigma,
    FALSE
  )$bound)
}

fit_from_beta <- function(case, beta) {
  control <- list(
    custom.beta = list(log(beta)),
    incremental.seed = stochastic_seed,
    incremental.final.theta = TRUE
  )
  elapsed <- system.time(model <- stm(
    case$documents, case$vocab, K = case$K,
    init.type = "Custom", inference = "Incremental",
    seed = model_seed, max.em.its = 100L, emtol = 1e-5,
    control = control, verbose = FALSE
  ))[["elapsed"]]
  list(
    beta = exp(model$beta$logbeta[[1L]]),
    bound = fresh_bound(model, case$documents),
    elapsed = unname(elapsed),
    epochs = model$convergence$its,
    converged = isTRUE(model$convergence$converged)
  )
}

# Self-contained shortest augmenting-path solver for one-to-one topic matching.
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

topic_metrics <- function(reference, candidate, top_n = 20L) {
  affinity <- sqrt(reference) %*% t(sqrt(candidate))
  assignment <- solve_assignment(-affinity)
  matched <- candidate[assignment, , drop = FALSE]
  matched_affinity <- rowSums(sqrt(reference * matched))
  top_overlap <- vapply(seq_len(nrow(reference)), function(k) {
    reference_top <- order(reference[k, ], decreasing = TRUE)[seq_len(top_n)]
    candidate_top <- order(matched[k, ], decreasing = TRUE)[seq_len(top_n)]
    length(intersect(reference_top, candidate_top))
  }, numeric(1L))
  c(
    mean_hellinger = mean(sqrt(pmax(0, 1 - matched_affinity))),
    max_hellinger = max(sqrt(pmax(0, 1 - matched_affinity))),
    mean_top20_overlap = mean(top_overlap)
  )
}

pairwise_stability <- function(fits, k) {
  pairs <- utils::combn(seq_along(fits), 2L)
  metrics <- apply(pairs, 2L, function(pair) {
    initial <- topic_metrics(
      fits[[pair[[1L]]]]$initial_beta,
      fits[[pair[[2L]]]]$initial_beta
    )
    final <- topic_metrics(
      fits[[pair[[1L]]]]$final_beta,
      fits[[pair[[2L]]]]$final_beta
    )
    c(
      initial_mean_hellinger = initial[["mean_hellinger"]],
      final_mean_hellinger = final[["mean_hellinger"]],
      initial_top20_overlap = initial[["mean_top20_overlap"]],
      final_top20_overlap = final[["mean_top20_overlap"]],
      anchor_overlap = length(intersect(
        fits[[pair[[1L]]]]$anchors,
        fits[[pair[[2L]]]]$anchors
      )) / k
    )
  })
  rowMeans(metrics)
}

raw_rows <- list()
stability_rows <- list()
fit_store <- list()
exact_store <- list()

for (case_name in names(cases)) {
  case <- cases[[case_name]]
  message(sprintf("Computing exact reference for %s", case_name))
  exact_initial <- exact_initialization(case)
  exact_final <- fit_from_beta(case, exact_initial$beta)
  exact_store[[case_name]] <- list(initial = exact_initial, final = exact_final)

  case_configs <- configurations
  if (mode == "full" && length(case$documents) >= 5000L) {
    case_configs <- case_configs[case_configs$large_case, , drop = FALSE]
  }
  fit_store[[case_name]] <- list()

  for (config_index in seq_len(nrow(case_configs))) {
    config <- case_configs[config_index, , drop = FALSE]
    config_fits <- vector("list", length(rp_seeds))
    for (seed_index in seq_along(rp_seeds)) {
      rp_seed <- rp_seeds[[seed_index]]
      message(sprintf(
        "Fitting %-14s %-11s rp.seed=%-4d p=%-4d density=%.2f",
        case_name, config$config, rp_seed, config$p, config$density
      ))
      initialized <- projected_initialization(
        case, config$p, config$density, rp_seed
      )
      final <- fit_from_beta(case, initialized$beta)
      initial_exact <- topic_metrics(exact_initial$beta, initialized$beta)
      final_exact <- topic_metrics(exact_final$beta, final$beta)
      anchor_exact_overlap <- length(intersect(
        exact_initial$anchors, initialized$anchors
      )) / case$K

      config_fits[[seed_index]] <- list(
        initial_beta = initialized$beta,
        final_beta = final$beta,
        anchors = initialized$anchors
      )
      raw_rows[[length(raw_rows) + 1L]] <- data.frame(
        case = case_name,
        config = config$config,
        p = config$p,
        density = config$density,
        nonzeros = config$nonzeros,
        rp_seed = rp_seed,
        init_seconds = initialized$elapsed,
        final_seconds = final$elapsed,
        final_epochs = final$epochs,
        final_converged = final$converged,
        final_bound = final$bound,
        bound_per_token = final$bound / case$ntokens,
        bound_delta_exact = final$bound - exact_final$bound,
        exact_anchor_overlap = anchor_exact_overlap,
        initial_exact_hellinger = initial_exact[["mean_hellinger"]],
        final_exact_hellinger = final_exact[["mean_hellinger"]],
        initial_exact_top20 = initial_exact[["mean_top20_overlap"]],
        final_exact_top20 = final_exact[["mean_top20_overlap"]]
      )
    }
    fit_store[[case_name]][[config$config]] <- config_fits
    stability <- pairwise_stability(config_fits, case$K)
    stability_rows[[length(stability_rows) + 1L]] <- data.frame(
      case = case_name,
      config = config$config,
      initial_pairwise_hellinger =
        stability[["initial_mean_hellinger"]],
      final_pairwise_hellinger = stability[["final_mean_hellinger"]],
      initial_pairwise_top20 = stability[["initial_top20_overlap"]],
      final_pairwise_top20 = stability[["final_top20_overlap"]],
      pairwise_anchor_overlap = stability[["anchor_overlap"]]
    )
  }
}

raw <- do.call(rbind, raw_rows)
stability <- do.call(rbind, stability_rows)
rownames(raw) <- rownames(stability) <- NULL

baseline_raw <- if (mode == "full") raw else readRDS(baseline_file)$raw
default_bound <- setNames(vapply(names(cases), function(case_name) {
  baseline_raw$final_bound[
    baseline_raw$case == case_name & baseline_raw$config == "p256-s05" &
      baseline_raw$rp_seed == 9L
  ][[1L]]
}, numeric(1L)), names(cases))
raw$bound_delta_default <- raw$final_bound - default_bound[raw$case]

summary <- do.call(rbind, lapply(split(raw, interaction(raw$case, raw$config,
                                                        drop = TRUE)),
                                  function(x) {
  data.frame(
    case = x$case[[1L]],
    config = x$config[[1L]],
    p = x$p[[1L]],
    density = x$density[[1L]],
    nonzeros = x$nonzeros[[1L]],
    nseeds = nrow(x),
    converged = sum(x$final_converged),
    median_init_seconds = median(x$init_seconds),
    median_final_seconds = median(x$final_seconds),
    median_bound = median(x$final_bound),
    bound_sd = sd(x$final_bound),
    bound_range_per_token = diff(range(x$bound_per_token)),
    median_delta_exact = median(x$bound_delta_exact),
    best_delta_exact = max(x$bound_delta_exact),
    worst_delta_exact = min(x$bound_delta_exact),
    beats_exact = sum(x$bound_delta_exact >= 0),
    median_delta_default = median(x$bound_delta_default),
    best_delta_default = max(x$bound_delta_default),
    fixed_seed9_delta_exact = x$bound_delta_exact[x$rp_seed == 9L][[1L]],
    mean_exact_anchor_overlap = mean(x$exact_anchor_overlap),
    mean_initial_exact_hellinger = mean(x$initial_exact_hellinger),
    mean_final_exact_hellinger = mean(x$final_exact_hellinger)
  )
}))
rownames(summary) <- NULL
summary <- merge(summary, stability, by = c("case", "config"), all.x = TRUE)
summary <- summary[order(match(summary$case, names(cases)),
                         match(summary$config, configurations$config)), ]
rownames(summary) <- NULL

exact_summary <- do.call(rbind, lapply(names(cases), function(case_name) {
  case <- cases[[case_name]]
  exact <- exact_store[[case_name]]
  data.frame(
    case = case_name,
    init_seconds = exact$initial$elapsed,
    final_seconds = exact$final$elapsed,
    final_epochs = exact$final$epochs,
    final_converged = exact$final$converged,
    final_bound = exact$final$bound,
    bound_per_token = exact$final$bound / case$ntokens
  )
}))
rownames(exact_summary) <- NULL

plot_file <- sub("\\.rds$", ".png", output, ignore.case = TRUE)
grDevices::png(plot_file, width = 1500L, height = 1100L, res = 180L)
old_par <- graphics::par(mfrow = c(length(cases), 2L),
                         mar = c(4.5, 4.5, 3, 1))
for (case_name in names(cases)) {
  x <- summary[summary$case == case_name, , drop = FALSE]
  graphics::plot(
    x$initial_pairwise_hellinger, x$median_delta_exact,
    pch = 19L, col = "#2166AC",
    xlab = "pairwise initialization Hellinger",
    ylab = "median final-bound delta vs exact",
    main = case_name
  )
  graphics::text(x$initial_pairwise_hellinger, x$median_delta_exact,
                 labels = x$config, pos = 3L, cex = .65)
  graphics::abline(h = 0, col = "#999999", lty = 3L)

  graphics::plot(
    x$final_pairwise_hellinger, x$bound_range_per_token,
    pch = 19L, col = "#B2182B",
    xlab = "pairwise final-topic Hellinger",
    ylab = "final-bound range per token",
    main = case_name
  )
  graphics::text(x$final_pairwise_hellinger, x$bound_range_per_token,
                 labels = x$config, pos = 3L, cex = .65)
}
graphics::par(old_par)
grDevices::dev.off()

result <- list(
  generated = Sys.time(),
  stm_version = as.character(packageVersion("stm")),
  model_seed = model_seed,
  stochastic_seed = stochastic_seed,
  rp_seeds = rp_seeds,
  configurations = configurations,
  exact = exact_summary,
  raw = raw,
  summary = summary,
  plot_file = normalizePath(plot_file, mustWork = FALSE)
)
saveRDS(result, output)
utils::write.csv(raw, sub("\\.rds$", "-raw.csv", output), row.names = FALSE)
utils::write.csv(summary, sub("\\.rds$", "-summary.csv", output),
                 row.names = FALSE)

cat("\nExact references:\n")
print(exact_summary, row.names = FALSE, digits = 6)
cat("\nSketch grid summary:\n")
print(summary, row.names = FALSE, digits = 6)
cat("\nSaved results to", normalizePath(output, mustWork = FALSE), "\n")
cat("Saved diagnostic plot to", normalizePath(plot_file, mustWork = FALSE), "\n")
