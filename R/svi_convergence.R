# Convergence Checking for Stochastic Variational Inference
#
# Mini-batch ELBO is noisy, so we use windowed averaging and
# patience-based early stopping rather than strict tolerance checking.

#' Initialize SVI Convergence Tracking
#'
#' Creates the convergence structure for tracking ELBO during
#' stochastic variational inference.
#'
#' @param settings Settings list containing max_epochs and patience
#'
#' @return Convergence tracking list
#'
#' @keywords internal
initialize_svi_convergence <- function(settings) {
  list(
    bound_history = numeric(0),      # Per-token mini-batch bounds (scaled)
    bound_window = numeric(0),       # Windowed per-token averages
    holdout_bound_history = numeric(0), # Holdout per-token ELBO history
    holdout_window = numeric(0),     # Holdout windowed per-token averages
    best_window = -Inf,              # Best windowed ELBO seen
    patience_counter = 0,            # Iterations without improvement
    train_stop = FALSE,              # Should stop based on training?
    holdout_best = -Inf,             # Best holdout windowed ELBO
    holdout_patience_counter = 0,    # Holdout patience counter
    holdout_stop = FALSE,            # Should stop based on holdout?
    its = 0,                         # Total iterations
    epoch = 0,                       # Current epoch (approximate)
    converged = FALSE,               # Did we converge?
    stopits = FALSE                  # Should we stop iterating?
  )
}

#' Check SVI Convergence
#'
#' Checks whether stochastic VI has converged using windowed ELBO
#' averaging and patience-based early stopping. Mini-batch ELBO is
#' noisy, so we smooth over a window of recent iterations.
#'
#' @param batch_elbo ELBO from current mini-batch (sum of suffstats$bound)
#' @param convergence Convergence tracking structure
#' @param batch_size Current mini-batch size
#' @param N Total number of documents
#' @param settings Settings list with convergence parameters
#'
#' @return Updated convergence structure
#'
#' @keywords internal
svi_convergence_check <- function(batch_elbo, convergence,
                                  batch_size, N, settings) {
  # Scale bound to full-data equivalent and normalize by tokens
  # This makes bounds comparable across different batch sizes
  scaled_bound <- batch_elbo * (N / batch_size)
  train_tokens <- if(!is.null(settings$svi$train_tokens)) {
    settings$svi$train_tokens
  } else if(!is.null(settings$dim$wcounts$x)) {
    sum(settings$dim$wcounts$x)
  } else {
    N
  }
  scaled_bound <- scaled_bound / train_tokens

  if(!is.finite(scaled_bound)) {
    warning("Non-finite ELBO encountered; skipping convergence update for this iteration.")
    scaled_bound <- NA_real_
  }

  # Add to history
  convergence$bound_history <- c(convergence$bound_history, scaled_bound)
  convergence$its <- convergence$its + 1

  # Update epoch counter (approximate - based on cumulative documents seen)
  convergence$epoch <- convergence$its * batch_size / N

  # Compute windowed average (last 10 iterations by default)
  window_size <- if(!is.null(settings$svi$convergence_window)) {
    settings$svi$convergence_window
  } else {
    10
  }

  if(length(convergence$bound_history) >= window_size) {
    recent <- tail(convergence$bound_history, window_size)
    if(all(!is.finite(recent))) {
      return(convergence)
    }
    window_mean <- mean(recent, na.rm=TRUE)
    window_var <- var(recent, na.rm=TRUE)
    if(!is.finite(window_mean)) {
      return(convergence)
    }

    convergence$bound_window <- c(convergence$bound_window, window_mean)

    # Check for improvement
    # We consider it an improvement if window_mean increases by at least
    # a small amount (to avoid stopping on numerical noise)
    improvement_tol <- if(!is.null(settings$svi$improvement_tol)) {
      settings$svi$improvement_tol
    } else {
      1e-4
    }
    if(!is.finite(improvement_tol)) {
      improvement_tol <- 1e-4
    }
    if(!is.finite(convergence$best_window)) {
      convergence$best_window <- -Inf
    }
    improvement_threshold <- improvement_tol * max(1, abs(convergence$best_window))
    if(!is.finite(improvement_threshold)) {
      improvement_threshold <- 0
    }
    if(!is.finite(convergence$best_window)) {
      improved <- is.finite(window_mean)
    } else {
      improved <- isTRUE(
        is.finite(window_mean) && is.finite(improvement_threshold) &&
        window_mean > convergence$best_window + improvement_threshold
      )
    }
    if(improved) {
      convergence$best_window <- window_mean
      convergence$patience_counter <- 0
    } else {
      convergence$patience_counter <- convergence$patience_counter + 1
    }

    # Early stopping via patience
    patience <- if(!is.null(settings$svi$patience)) {
      settings$svi$patience
    } else {
      20
    }
    if(!is.finite(patience)) {
      patience <- 20
    }

    if(convergence$patience_counter >= patience) {
      convergence$train_stop <- TRUE
      early_stop <- if(!is.null(settings$svi$early_stop)) {
        settings$svi$early_stop
      } else {
        "train"
      }
      if(early_stop == "train") {
        convergence$converged <- TRUE
        convergence$stopits <- TRUE
        if(settings$verbose) {
          cat(sprintf("\nConverged: no improvement for %d iterations (%.1f epochs)\n",
                     patience, convergence$epoch))
          cat(sprintf("Final windowed ELBO (per token): %.4f\n", window_mean))
        }
      }
    }
  }

  # Max iterations check
  max_iters <- settings$svi$max_iters
  if(!is.null(max_iters) && convergence$its >= max_iters) {
    convergence$stopits <- TRUE
    if(settings$verbose) {
      cat(sprintf("\nReached maximum iterations: %d\n", max_iters))
    }
  }

  # Max epochs check
  max_epochs <- settings$svi$max_epochs
  if(!is.null(max_epochs) && convergence$epoch >= max_epochs) {
    convergence$stopits <- TRUE
    if(settings$verbose) {
      cat(sprintf("\nReached maximum epochs: %.1f\n", max_epochs))
    }
  }

  return(convergence)
}

#' Report SVI Progress
#'
#' Reports iteration progress during stochastic variational inference.
#' Similar to report() in STMreport.R but adapted for mini-batch updates.
#'
#' @param iter Current iteration number
#' @param convergence Convergence tracking structure
#' @param settings Settings list
#'
#' @keywords internal
report_svi_progress <- function(iter, convergence, settings) {
  verbose <- settings$verbose
  reportevery <- settings$topicreportevery

  # Only report at specified intervals
  if(!verbose || (iter != 1 && iter %% reportevery != 0)) {
    return(invisible())
  }

  # Get current bounds
  recent_elbo <- convergence$bound_history[length(convergence$bound_history)]
  windowed_elbo <- ifelse(length(convergence$bound_window) > 0,
                          tail(convergence$bound_window, 1), NA)

  # During warm-up (first few iterations before window fills)
  if(iter <= settings$svi$convergence_window) {
    msg <- sprintf("Completing Iteration %d (approx. per word bound = %.3f) \n",
                   iter, recent_elbo)
  } else {
    # After warm-up: show windowed bound and convergence info
    patience <- convergence$patience_counter
    max_patience <- ifelse(is.null(settings$svi$patience), 20, settings$svi$patience)

    msg <- sprintf("Completing Iteration %d (approx. per word bound = %.3f, patience = %d/%d) \n",
                   iter, windowed_elbo, patience, max_patience)
  }

  cat(msg)

  # Topic preview at specified intervals (match base stm behavior)
  # Only show topics periodically, not every iteration
  if(iter > 1 && iter %% (reportevery * 5) == 0) {
    # Similar to STMreport.R logic
    beta <- settings$beta_current  # Need to pass this in
    vocab <- settings$vocab

    if(!is.null(beta) && !is.null(vocab)) {
      wordmat <- apply(beta[[1]], 1, function(x) vocab[order(x, decreasing=TRUE)[1:5]])
      labs <- apply(wordmat, 2, function(x) paste(x, collapse=", "))
      toprint <- sprintf("     Topic %i: %s \n", 1:length(labs), labs)
      cat(toprint)
    }
  }
}

#' Compute Full ELBO for Validation
#'
#' Computes ELBO on the entire dataset (expensive!). Useful for
#' periodic validation during training to see true convergence.
#'
#' @param documents Full document list
#' @param mu Mu parameter list
#' @param sigma Covariance matrix
#' @param beta Beta parameter list
#' @param lambda Current lambda values
#' @param betaindex Content covariate indices
#' @param settings Settings list
#'
#' @return Full-dataset ELBO
#'
#' @keywords internal
compute_full_elbo <- function(documents, mu, sigma, beta, lambda,
                              betaindex, settings) {
  # Run E-step on full dataset (no updates, just compute bound)
  update_mu <- !is.null(settings$prevalence)
  estep_fn <- if(exists("estep", mode="function")) estep else stm:::estep
  suffstats <- estep_fn(
    documents = documents,
    beta.index = betaindex,
    update.mu = update_mu,
    beta = beta,
    lambda.old = lambda,
    mu = mu,
    sigma = sigma,
    verbose = FALSE
  )

  return(sum(suffstats$bound))
}
