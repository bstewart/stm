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
    bound_history = numeric(0),      # All mini-batch bounds (scaled)
    bound_window = numeric(0),       # Windowed averages
    best_window = -Inf,              # Best windowed ELBO seen
    patience_counter = 0,            # Iterations without improvement
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
  # Scale bound to full-data equivalent
  # This makes bounds comparable across different batch sizes
  scaled_bound <- batch_elbo * (N / batch_size)

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
    window_mean <- mean(recent)
    window_var <- var(recent)

    convergence$bound_window <- c(convergence$bound_window, window_mean)

    # Check for improvement
    # We consider it an improvement if window_mean increases by at least
    # a small amount (to avoid stopping on numerical noise)
    improvement_threshold <- 1e-3
    if(window_mean > convergence$best_window + improvement_threshold) {
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

    if(convergence$patience_counter >= patience) {
      convergence$converged <- TRUE
      convergence$stopits <- TRUE
      if(settings$verbose) {
        cat(sprintf("\nConverged: no improvement for %d iterations (%.1f epochs)\n",
                   patience, convergence$epoch))
        cat(sprintf("Final windowed ELBO: %.2f\n", window_mean))
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

  # Compute per-word bound (comparable to standard stm)
  ntokens <- sum(settings$dim$wcounts$x)
  tokenll <- recent_elbo / ntokens

  # During warm-up (first few iterations before window fills)
  if(iter <= settings$svi$convergence_window) {
    msg <- sprintf("Completing Iteration %d (approx. per word bound = %.3f) \n",
                   iter, tokenll)
  } else {
    # After warm-up: show windowed bound and convergence info
    window_tokenll <- windowed_elbo / ntokens
    patience <- convergence$patience_counter
    max_patience <- ifelse(is.null(settings$svi$patience), 20, settings$svi$patience)

    msg <- sprintf("Completing Iteration %d (approx. per word bound = %.3f, patience = %d/%d) \n",
                   iter, window_tokenll, patience, max_patience)
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
  suffstats <- estep(
    documents = documents,
    beta.index = betaindex,
    update.mu = FALSE,
    beta = beta,
    lambda.old = lambda,
    mu = mu,
    sigma = sigma,
    verbose = FALSE
  )

  return(sum(suffstats$bound))
}
