# Gradient Computation for Stochastic Variational Inference
#
# This file converts sufficient statistics from the E-step into
#scaled gradients for stochastic optimization. The critical scaling
# factor N/batch_size ensures unbiased gradient estimates.

#' Compute SVI Gradients from Sufficient Statistics
#'
#' Converts E-step sufficient statistics to properly scaled gradients
#' for mini-batch stochastic optimization. The scaling factor N/batch_size
#' is CRITICAL for correctness - it converts mini-batch statistics into
#' unbiased estimates of full-dataset gradients.
#'
#' @param suffstats Sufficient statistics from E-step (lambda, sigma.ss, beta.ss, bound)
#' @param N Total number of documents in corpus
#' @param batch_size Number of documents in current mini-batch
#' @param mu List with $mu matrix of prevalence parameters
#' @param sigma Covariance matrix
#' @param beta List with $beta list of topic-word distributions
#' @param settings Settings list from stm
#'
#' @return List of gradients (mu, sigma, log_beta)
#'
#' @keywords internal
compute_svi_gradients <- function(suffstats, N, batch_size,
                                  mu, sigma, beta, settings) {
  # Critical scaling factor: converts mini-batch stats to full-data gradient estimates
  scale_factor <- N / batch_size

  K <- nrow(beta[[1]])
  Kminus1 <- K - 1

  # --- Gradient for mu (prevalence parameters) ---
  # NOTE: When prevalence covariates are present (mu is N×K-1), we skip gradient computation
  # because mu = X %*% gamma and gamma is updated separately via opt.mu()
  # Only compute gradients for CTM case (no covariates, mu is scalar/global)
  if(!is.null(settings$prevalence)) {
    # Skip mu gradient - will be updated via gamma
    grad_mu <- NULL
  } else {
    # For CTM case (no covariates): gradient is scaled mean difference
    # The sufficient statistic is lambda (document-level variational parameters)
    # Gradient of ELBO w.r.t. mu: sum_d (lambda_d - mu)
    grad_mu <- scale_factor * (colMeans(suffstats$lambda) - as.numeric(mu))
    grad_mu <- matrix(grad_mu, nrow=nrow(mu), ncol=ncol(mu))
  }

  # --- Gradient for sigma (covariance matrix) ---
  # Center lambda around mu
  # When mu is document-specific (prevalence covariates), we can't use simple centering
  # Instead, skip centering for now (will be handled during periodic gamma updates)
  if(!is.null(settings$prevalence)) {
    # For document-specific mu, don't center (approximation)
    lambda_centered <- suffstats$lambda
  } else {
    # For global mu, center normally
    lambda_centered <- sweep(suffstats$lambda, 2,
                            as.numeric(mu), FUN="-")
  }

  # Empirical covariance from batch
  empirical_cov <- crossprod(lambda_centered) / batch_size

  # Posterior variance from E-step (nu term)
  # suffstats$sigma is the sum of document-level posterior variances
  scaled_sigma_ss <- suffstats$sigma / batch_size

  # Natural gradient for covariance (uses inverse)
  # Gradient of ELBO w.r.t. sigma:
  # -0.5 * N * sigma^-1 + 0.5 * sum_n [(lambda_n - mu)(lambda_n - mu)' + nu_n] * sigma^-1
  siginv <- solve(sigma)
  grad_sigma <- scale_factor * 0.5 * siginv %*%
                (empirical_cov + scaled_sigma_ss) %*% siginv
  grad_sigma <- grad_sigma - 0.5 * scale_factor * siginv

  # Make symmetric (can lose symmetry due to numerical issues)
  grad_sigma <- (grad_sigma + t(grad_sigma)) / 2

  # Apply sigma prior (regularization toward diagonal)
  if(settings$sigma$prior > 0) {
    sigma_diag <- diag(diag(sigma))
    # Prior gradient pulls sigma toward diagonal
    grad_sigma <- grad_sigma -
                  settings$sigma$prior * scale_factor * (sigma - sigma_diag)
  }

  # --- Gradient for beta (topic-word distributions) ---
  # Sufficient statistic beta.ss is the expected word counts under phi
  # For standard LDA-style beta:
  #   Gradient of ELBO w.r.t. log(beta_kv) = E[n_kv] - beta_kv * sum_v E[n_kv]
  # where E[n_kv] is the expected count of word v in topic k

  grad_log_beta <- vector("list", length(beta))
  for(a in 1:length(beta)) {
    # Scale sufficient statistics to estimate full-data counts
    scaled_counts <- scale_factor * suffstats$beta[[a]]

    # For each topic k, the gradient has two terms:
    # 1. Expected counts (positive - encourages high beta for frequent words)
    # 2. Normalization term (negative - enforces simplex constraint)

    # Gradient w.r.t. log(beta): E[counts] - beta * sum_v E[counts]
    # The second term enforces that sum_v beta_kv = 1
    grad_log_beta[[a]] <- scaled_counts -
                          beta[[a]] * rowSums(scaled_counts)
  }

  return(list(mu=grad_mu, sigma=grad_sigma, log_beta=grad_log_beta))
}

#' Clip Gradients to Prevent Divergence
#'
#' Applies element-wise clipping to gradients to prevent large
#' updates that can cause numerical instability or divergence.
#'
#' @param gradients List of gradients (mu, sigma, log_beta)
#' @param clip_value Maximum absolute value for gradient elements
#'
#' @return Clipped gradients
#'
#' @keywords internal
clip_gradients <- function(gradients, clip_value=5.0) {
  gradients$mu <- pmin(pmax(gradients$mu, -clip_value), clip_value)
  gradients$sigma <- pmin(pmax(gradients$sigma, -clip_value), clip_value)

  for(a in 1:length(gradients$log_beta)) {
    gradients$log_beta[[a]] <- pmin(pmax(gradients$log_beta[[a]],
                                         -clip_value), clip_value)
  }

  return(gradients)
}

#' Check for NaN or Inf in Gradients
#'
#' Detects numerical issues in computed gradients. If found,
#' prints diagnostic information and optionally stops execution.
#'
#' @param gradients List of gradients to check
#' @param iter Current iteration number (for error messages)
#' @param stop_on_error If TRUE, stops execution when NaN/Inf found
#'
#' @return TRUE if gradients are valid, FALSE otherwise
#'
#' @keywords internal
check_gradients <- function(gradients, iter, stop_on_error=TRUE) {
  has_nan <- FALSE
  msg <- ""

  # Check mu
  if(any(is.na(gradients$mu)) || any(is.infinite(gradients$mu))) {
    has_nan <- TRUE
    msg <- paste0(msg, "  Gradient for mu contains NaN or Inf\n")
  }

  # Check sigma
  if(any(is.na(gradients$sigma)) || any(is.infinite(gradients$sigma))) {
    has_nan <- TRUE
    msg <- paste0(msg, "  Gradient for sigma contains NaN or Inf\n")
  }

  # Check beta
  for(a in 1:length(gradients$log_beta)) {
    if(any(is.na(gradients$log_beta[[a]])) ||
       any(is.infinite(gradients$log_beta[[a]]))) {
      has_nan <- TRUE
      msg <- paste0(msg, sprintf("  Gradient for beta (aspect %d) contains NaN or Inf\n", a))
    }
  }

  if(has_nan) {
    full_msg <- sprintf("Iteration %d: Invalid gradients detected:\n%s", iter, msg)
    if(stop_on_error) {
      stop(full_msg)
    } else {
      warning(full_msg)
    }
    return(FALSE)
  }

  return(TRUE)
}
