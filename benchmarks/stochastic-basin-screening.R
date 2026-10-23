#!/usr/bin/env Rscript

# Can early variational bounds identify which randomized spectral sketch will
# converge to the best stochastic-inference basin?
#
# For each rp.seed this benchmark records:
#   initial:  the exact full-corpus E-step bound at spectral initialization;
#   one_step: a fresh full-corpus bound after the first exact E/M step; and
#   final:    a fresh full-corpus bound after stochastic convergence.
#
# The stochastic update seed is held fixed, so variation is attributable to
# the spectral sketch.  Epoch 1 of stochastic inference populates the complete
# sufficient-statistic cache before applying one batch-equivalent M-step;
# shuffled block replacement does not begin until epoch 2.
#
# Usage:
#   Rscript benchmarks/stochastic-basin-screening.R [nstarts] [output.rds]

suppressPackageStartupMessages(library(stm))

args <- commandArgs(trailingOnly = TRUE)
nstarts <- if (length(args)) as.integer(args[[1L]]) else 12L
if (is.na(nstarts) || nstarts < 5L || nstarts > 20L) {
  stop("nstarts must be between 5 and 20")
}
output <- if (length(args) > 1L) args[[2L]] else
  file.path(tempdir(), "stm-stochastic-basin-screening.rds")

available_rp_seeds <- c(
  9L, 17L, 41L, 89L, 151L, 3L, 5L, 7L, 11L, 13L,
  19L, 23L, 29L, 31L, 37L, 43L, 47L, 53L, 59L, 61L
)
rp_seeds <- available_rp_seeds[seq_len(nstarts)]
model_seed <- 91L
stochastic_seed <- model_seed + 104729L

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

fit_case <- function(case, rp_seed, one_step = FALSE) {
  control <- list(
    rp.seed = as.integer(rp_seed),
    incremental.seed = stochastic_seed,
    incremental.final.theta = TRUE
  )
  if (one_step) control$incremental.max.epochs <- 1L
  fit_args <- list(
    documents = case$documents,
    vocab = case$vocab,
    K = case$K,
    data = case$data,
    init.type = "Spectral",
    inference = "Incremental",
    seed = model_seed,
    max.em.its = 75L,
    emtol = 1e-5,
    control = control,
    verbose = FALSE
  )
  if (!is.null(case$prevalence)) fit_args$prevalence <- case$prevalence
  elapsed <- system.time(model <- do.call(stm, fit_args))[["elapsed"]]
  list(
    model = model,
    elapsed = unname(elapsed),
    fresh_bound = fresh_bound(model, case$documents)
  )
}

run_start <- function(case, rp_seed) {
  message(sprintf("Fitting %-22s rp.seed=%-4d checkpoint=final",
                  case$name, rp_seed))
  final <- fit_case(case, rp_seed, one_step = FALSE)
  message(sprintf("Fitting %-22s rp.seed=%-4d checkpoint=one-step",
                  case$name, rp_seed))
  one_step <- fit_case(case, rp_seed, one_step = TRUE)

  initial_final_run <- final$model$convergence$bound[[1L]]
  initial_one_step_run <- one_step$model$convergence$bound[[1L]]
  data.frame(
    case = case$name,
    rp_seed = rp_seed,
    ntokens = case$ntokens,
    initial_bound = initial_final_run,
    initial_replicate_difference =
      initial_one_step_run - initial_final_run,
    one_step_bound = one_step$fresh_bound,
    final_bound = final$fresh_bound,
    final_epochs = final$model$convergence$its,
    final_converged = isTRUE(final$model$convergence$converged),
    one_step_seconds = one_step$elapsed,
    final_seconds = final$elapsed
  )
}

raw <- do.call(rbind, lapply(cases, function(case) {
  do.call(rbind, lapply(rp_seeds, function(seed) run_start(case, seed)))
}))
rownames(raw) <- NULL

screening_summary <- function(x, checkpoint) {
  early <- x[[paste0(checkpoint, "_bound")]]
  final <- x$final_bound
  selected <- which.max(early)
  final_order <- order(final, decreasing = TRUE)
  early_order <- order(early, decreasing = TRUE)
  selected_rank <- match(selected, final_order)
  top_n <- min(3L, nrow(x))
  early_top <- early_order[seq_len(top_n)]
  final_top <- final_order[seq_len(top_n)]
  best_final_index <- final_order[[1L]]
  best_advanced_index <- early_top[[which.max(final[early_top])]]
  pairs <- utils::combn(seq_len(nrow(x)), 2L)
  early_direction <- sign(early[pairs[1L, ]] - early[pairs[2L, ]])
  final_direction <- sign(final[pairs[1L, ]] - final[pairs[2L, ]])
  usable <- early_direction != 0 & final_direction != 0
  spearman <- suppressWarnings(cor.test(early, final, method = "spearman",
                                        exact = FALSE))
  data.frame(
    case = x$case[[1L]],
    checkpoint = checkpoint,
    nstarts = nrow(x),
    converged = sum(x$final_converged),
    pearson = cor(early, final, method = "pearson"),
    spearman = unname(spearman$estimate),
    spearman_p = spearman$p.value,
    pairwise_order_accuracy = mean(early_direction[usable] ==
                                     final_direction[usable]),
    selected_seed = x$rp_seed[[selected]],
    selected_final_rank = selected_rank,
    best_final_seed = x$rp_seed[[best_final_index]],
    best_final_early_rank = match(best_final_index, early_order),
    selection_regret = max(final) - final[[selected]],
    selection_regret_per_token =
      (max(final) - final[[selected]]) / x$ntokens[[1L]],
    top3_overlap = length(intersect(early_top, final_top)),
    top3_recall = length(intersect(early_top, final_top)) / top_n,
    best_final_in_early_top3 = best_final_index %in% early_top,
    top3_screen_regret = max(final) - final[[best_advanced_index]],
    top3_screen_regret_per_token =
      (max(final) - final[[best_advanced_index]]) / x$ntokens[[1L]]
  )
}

summary <- do.call(rbind, lapply(split(raw, raw$case), function(x) {
  rbind(screening_summary(x, "initial"),
        screening_summary(x, "one_step"))
}))
rownames(summary) <- NULL

initial_reproducibility <- aggregate(
  abs(initial_replicate_difference) ~ case,
  data = raw,
  FUN = max
)
names(initial_reproducibility)[[2L]] <-
  "max_initial_bound_replicate_difference"

# A compact six-panel diagnostic.  Axes are per-token regret from the best
# value at each checkpoint, so zero is the upper-right corner in every panel.
plot_file <- if (grepl("\\.rds$", output, ignore.case = TRUE)) {
  sub("\\.rds$", ".png", output, ignore.case = TRUE)
} else {
  paste0(output, ".png")
}
grDevices::png(plot_file, width = 1500L, height = 1800L, res = 180L)
old_par <- graphics::par(mfrow = c(3L, 2L), mar = c(4.2, 4.4, 2.6, 1.0),
                         oma = c(0, 0, 2, 0))
for (case_name in names(cases)) {
  x <- raw[raw$case == case_name, , drop = FALSE]
  for (checkpoint in c("initial", "one_step")) {
    early <- x[[paste0(checkpoint, "_bound")]] / x$ntokens
    final <- x$final_bound / x$ntokens
    early <- early - max(early)
    final <- final - max(final)
    graphics::plot(
      early, final,
      pch = ifelse(x$rp_seed == 9L, 19L, 1L),
      col = ifelse(x$rp_seed == 9L, "#B2182B", "#2166AC"),
      xlab = sprintf("%s bound regret per token", checkpoint),
      ylab = "final bound regret per token",
      main = case_name
    )
    graphics::text(early, final, labels = x$rp_seed, pos = 3L,
                   cex = 0.55, col = "#333333")
    graphics::abline(h = 0, v = 0, col = "#BBBBBB", lty = 3L)
    fit <- stats::lm(final ~ early)
    graphics::abline(fit, col = "#666666", lty = 2L)
    rho <- summary$spearman[
      summary$case == case_name & summary$checkpoint == checkpoint
    ]
    graphics::legend("bottomright", legend = sprintf("Spearman rho = %.2f", rho),
                     bty = "n", cex = 0.75)
  }
}
graphics::mtext("Can an early bound screen randomized spectral starts?",
                outer = TRUE, line = 0.3, cex = 1.2)
graphics::par(old_par)
grDevices::dev.off()

result <- list(
  generated = Sys.time(),
  stm_version = as.character(packageVersion("stm")),
  model_seed = model_seed,
  stochastic_seed = stochastic_seed,
  rp_seeds = rp_seeds,
  raw = raw,
  summary = summary,
  initial_reproducibility = initial_reproducibility,
  plot_file = normalizePath(plot_file, mustWork = FALSE)
)
saveRDS(result, output)
utils::write.csv(raw, sub("\\.rds$", "-raw.csv", output), row.names = FALSE)
utils::write.csv(summary, sub("\\.rds$", "-summary.csv", output),
                 row.names = FALSE)

cat("\nEarly-screening summary:\n")
print(summary, row.names = FALSE, digits = 5)
cat("\nInitial-bound repeatability check:\n")
print(initial_reproducibility, row.names = FALSE, digits = 5)
cat("\nSaved results to", normalizePath(output, mustWork = FALSE), "\n")
cat("Saved diagnostic plot to", normalizePath(plot_file, mustWork = FALSE), "\n")
