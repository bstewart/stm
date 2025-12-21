# Adam Optimizer Implementation for Stochastic Variational Inference
#
# This file implements the Adam optimizer (Kingma & Ba, 2014) for
# updating STM parameters in stochastic variational inference.

#' Adam Optimization Step
#'
#' Performs one step of Adam optimization on model parameters.
#'
#' @param adam_state List containing Adam optimizer state (first/second moments)
#' @param gradients List of gradients for mu, sigma, beta, and gamma parameters
#' @param lr Learning rate (step size)
#' @param iter Current iteration number (for bias correction)
#' @param beta1 Exponential decay rate for first moment (default: 0.9)
#' @param beta2 Exponential decay rate for second moment (default: 0.999)
#' @param eps Small constant for numerical stability (default: 1e-8)
#'
#' @return Updated adam_state with computed parameter updates
#'
#' @keywords internal
adam_update_step <- function(adam_state, gradients, lr, iter,
                             beta1=0.9, beta2=0.999, eps=1e-8) {
  # Increment iteration counter
  adam_state$t <- iter

  # --- Update mu parameters ---
  # Skip if NULL (prevalence covariates case - gamma updated separately)
  if(!is.null(gradients$mu)) {
    adam_state$m_mu <- beta1 * adam_state$m_mu + (1 - beta1) * gradients$mu
    adam_state$v_mu <- beta2 * adam_state$v_mu + (1 - beta2) * gradients$mu^2

    # Bias correction
    m_mu_hat <- adam_state$m_mu / (1 - beta1^iter)
    v_mu_hat <- adam_state$v_mu / (1 - beta2^iter)

    # Compute update (NOT applied yet, returned for caller to apply)
    adam_state$update_mu <- lr * m_mu_hat / (sqrt(v_mu_hat) + eps)
  } else {
    # No gradient - don't update
    adam_state$update_mu <- NULL
  }

  # --- Update gamma parameters (prevalence regression) ---
  if(!is.null(gradients$gamma)) {
    adam_state$m_gamma <- beta1 * adam_state$m_gamma +
                          (1 - beta1) * gradients$gamma
    adam_state$v_gamma <- beta2 * adam_state$v_gamma +
                          (1 - beta2) * gradients$gamma^2

    m_gamma_hat <- adam_state$m_gamma / (1 - beta1^iter)
    v_gamma_hat <- adam_state$v_gamma / (1 - beta2^iter)

    adam_state$update_gamma <- lr * m_gamma_hat / (sqrt(v_gamma_hat) + eps)
  } else {
    adam_state$update_gamma <- NULL
  }

  # --- Update sigma parameters ---
  adam_state$m_sigma <- beta1 * adam_state$m_sigma +
                        (1 - beta1) * gradients$sigma
  adam_state$v_sigma <- beta2 * adam_state$v_sigma +
                        (1 - beta2) * gradients$sigma^2

  m_sigma_hat <- adam_state$m_sigma / (1 - beta1^iter)
  v_sigma_hat <- adam_state$v_sigma / (1 - beta2^iter)

  adam_state$update_sigma <- lr * m_sigma_hat / (sqrt(v_sigma_hat) + eps)

  # --- Update beta parameters (list of matrices) ---
  adam_state$update_beta <- vector("list", length(gradients$log_beta))
  for(a in 1:length(gradients$log_beta)) {
    adam_state$m_beta[[a]] <- beta1 * adam_state$m_beta[[a]] +
                              (1 - beta1) * gradients$log_beta[[a]]
    adam_state$v_beta[[a]] <- beta2 * adam_state$v_beta[[a]] +
                              (1 - beta2) * gradients$log_beta[[a]]^2

    m_beta_hat <- adam_state$m_beta[[a]] / (1 - beta1^iter)
    v_beta_hat <- adam_state$v_beta[[a]] / (1 - beta2^iter)

    adam_state$update_beta[[a]] <- lr * m_beta_hat /
                                   (sqrt(v_beta_hat) + eps)
  }

  return(adam_state)
}

#' Update Beta Parameters from Adam
#'
#' Applies Adam updates to beta (topic-word distributions) while maintaining
#' probability constraints (non-negativity and simplex).
#'
#' @param beta_list List of K x V matrices of topic-word probabilities
#' @param update_list List of update matrices from Adam optimizer
#'
#' @return Updated beta_list with probabilities properly normalized
#'
#' @keywords internal
update_beta_from_adam <- function(beta_list, update_list) {
  # Beta is optimized in log space for unconstrained optimization
  # Apply update to log(beta), then re-normalize
  for(a in 1:length(beta_list)) {
    # Convert to log space
    log_beta <- log(beta_list[[a]] + 1e-10)  # Add small constant to avoid log(0)

    # Apply Adam update
    log_beta <- log_beta + update_list[[a]]

    # Convert back to probability space via softmax (row-wise)
    # This ensures non-negativity and simplex constraint
    log_beta_max <- apply(log_beta, 1, max)  # For numerical stability
    log_beta_shifted <- sweep(log_beta, 1, log_beta_max, FUN="-")
    exp_beta <- exp(log_beta_shifted)
    beta_list[[a]] <- exp_beta / rowSums(exp_beta)

    # Handle any numerical issues
    beta_list[[a]][is.na(beta_list[[a]])] <- 1/ncol(beta_list[[a]])
    beta_list[[a]][beta_list[[a]] < 1e-10] <- 1e-10
    beta_list[[a]] <- beta_list[[a]] / rowSums(beta_list[[a]])
  }
  return(beta_list)
}

#' Ensure Sigma is Positive Definite
#'
#' After Adam updates, sigma may lose positive definiteness due to
#' numerical issues. This function projects it back to PD cone.
#'
#' @param sigma Covariance matrix (may not be PD)
#'
#' @return Corrected positive definite sigma matrix
#'
#' @keywords internal
ensure_sigma_pd <- function(sigma) {
  # Force symmetry
  sigma <- (sigma + t(sigma)) / 2

  # Check if positive definite
  min_eig <- tryCatch({
    min(eigen(sigma, symmetric=TRUE, only.values=TRUE)$values)
  }, error=function(e) {
    -1  # If eigen fails, assume not PD
  })

  # If not PD, add small ridge to diagonal
  if(min_eig < 1e-6) {
    ridge <- max(1e-6 - min_eig, 1e-8)
    sigma <- sigma + ridge * diag(nrow(sigma))

    # Re-check symmetry
    sigma <- (sigma + t(sigma)) / 2
  }

  return(sigma)
}

#' Initialize Adam Optimizer State
#'
#' Creates the state structure for Adam optimizer, including
#' first and second moment estimates for all parameters.
#'
#' @param mu List with $mu matrix (prevalence parameters)
#' @param sigma Covariance matrix
#' @param beta List with $beta list of topic-word matrices
#' @param gamma Prevalence regression coefficients (optional)
#' @param has_prevalence If TRUE, skip mu initialization (gamma updated separately)
#'
#' @return adam_state list with initialized moments
#'
#' @keywords internal
initialize_adam_state <- function(mu, sigma, beta, gamma=NULL, has_prevalence=FALSE) {
  # When prevalence covariates are present, skip mu Adam state (gamma updated separately)
  if(has_prevalence) {
    m_mu <- NULL
    v_mu <- NULL
  } else {
    m_mu <- matrix(0, nrow=nrow(mu), ncol=ncol(mu))
    v_mu <- matrix(0, nrow=nrow(mu), ncol=ncol(mu))
  }
  if(has_prevalence && !is.null(gamma)) {
    m_gamma <- matrix(0, nrow=nrow(gamma), ncol=ncol(gamma))
    v_gamma <- matrix(0, nrow=nrow(gamma), ncol=ncol(gamma))
  } else {
    m_gamma <- NULL
    v_gamma <- NULL
  }

  list(
    # Mu parameters (prevalence) - NULL if has_prevalence=TRUE
    m_mu = m_mu,
    v_mu = v_mu,

    # Gamma parameters (prevalence regression)
    m_gamma = m_gamma,
    v_gamma = v_gamma,

    # Sigma parameters (covariance)
    m_sigma = matrix(0, nrow=nrow(sigma), ncol=ncol(sigma)),
    v_sigma = matrix(0, nrow=nrow(sigma), ncol=ncol(sigma)),

    # Beta parameters (topic-word distributions)
    m_beta = lapply(beta, function(b) matrix(0, nrow=nrow(b), ncol=ncol(b))),
    v_beta = lapply(beta, function(b) matrix(0, nrow=nrow(b), ncol=ncol(b))),

    # Iteration counter
    t = 0,

    # Placeholders for updates (will be filled by adam_update_step)
    update_mu = NULL,
    update_sigma = NULL,
    update_beta = NULL,
    update_gamma = NULL
  )
}
