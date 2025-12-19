# Final Theta Computation for SVI Models
#
# This file provides functionality to compute accurate document-topic
# proportions (theta) via a final complete E-step pass through all documents.

#' Check if Model Has Theta Computed
#'
#' Internal helper to check if an STM model has theta computed and
#' provide a helpful error message if not.
#'
#' @param model STM model object
#' @param function_name Name of calling function (for error message)
#' @keywords internal
check_theta_computed <- function(model, function_name) {
  if(is.null(model$theta)) {
    stop(
      sprintf("Function '%s' requires model$theta, but it is NULL.\n\n", function_name),
      "This model was fit with stm_svi(..., compute_final_theta=FALSE) ",
      "to skip the final E-step.\n\n",
      "To compute theta, run:\n",
      "  model <- update_svi_theta(model, documents)\n\n",
      "See ?update_svi_theta for details."
    )
  }
}

#' Perform Final E-Step for SVI Model
#'
#' Internal function that performs a complete E-step pass through all
#' documents using final learned parameters. Used both at the end of
#' stm_svi() training and by update_svi_theta() for post-hoc updates.
#'
#' @param documents Document list
#' @param betaindex Beta index for each document
#' @param mu Final mu matrix (K-1 x 1)
#' @param sigma Final sigma matrix (K-1 x K-1)
#' @param beta Final beta list (list of K x V matrices)
#' @param lambda_init Initial lambda for warm start (N x K-1)
#' @param parallel_chunks Number of parallel chunks (1=serial)
#' @param verbose Print progress
#'
#' @return List with lambda (N x K-1), theta (N x K), and bound
#'
#' @keywords internal
perform_final_estep <- function(documents, betaindex, mu, sigma, beta,
                                lambda_init, parallel_chunks=1, verbose=TRUE) {
  N <- length(documents)

  # Determine if mu is document-specific
  # If mu is (K-1) × N matrix, then update_mu=TRUE
  # If mu is (K-1) × 1 matrix or vector, then update_mu=FALSE
  update_mu <- (is.matrix(mu) && ncol(mu) == N)

  # Call parallel E-step on all documents
  final_suffstats <- tryCatch({
    estep_parallel(
      documents = documents,
      doc_indices = 1:N,
      betaindex = betaindex,
      beta = beta,
      lambda_old = lambda_init,
      mu = mu,
      sigma = sigma,
      parallel_chunks = parallel_chunks,
      verbose = verbose,
      update_mu = update_mu
    )
  }, error = function(e) {
    stop(sprintf("Error in final E-step: %s", conditionMessage(e)))
  })

  # Extract fresh lambda
  lambda <- final_suffstats$lambda

  # Compute theta using numerically stable approach
  # Add reference category (K-th topic with lambda=0)
  lambda_full <- cbind(lambda, 0)

  # Use matrixStats::rowLogSumExps for numerical stability
  # theta_ik = exp(lambda_ik) / sum_j exp(lambda_ij)
  # This is equivalent to row.lse() but accessed directly
  log_sum_exp <- matrixStats::rowLogSumExps(lambda_full)
  theta <- exp(lambda_full - log_sum_exp)

  # Ensure proper normalization (should already be normalized, but verify)
  theta <- theta / rowSums(theta)

  return(list(
    lambda = lambda,
    theta = theta,
    bound = final_suffstats$bound
  ))
}


#' Update Theta for SVI Model
#'
#' Computes document-topic proportions (theta) for an SVI model that was fit
#' without computing final theta. Performs a complete E-step on all documents
#' using the final parameter values.
#'
#' When \code{stm_svi()} is called with \code{compute_final_theta=FALSE}, the
#' returned model will have \code{theta=NULL}. This function computes theta
#' by performing a final E-step on all documents. This is useful for:
#' \itemize{
#'   \item Very large corpora where you want to avoid the final E-step during training
#'   \item Computing theta only when needed for specific analyses
#'   \item Updating theta after model parameters have been modified
#' }
#'
#' @param model An STM model object fitted with \code{stm_svi()}
#' @param documents Original document list (if not stored in model)
#' @param cores Number of CPU cores for parallel processing (NULL=auto-detect, 1=serial)
#' @param verbose Print progress messages
#'
#' @return Updated model with \code{theta} and \code{eta} computed. The
#'   \code{convergence$theta_computed} flag is set to TRUE.
#'
#' @details
#' The function performs a complete E-step pass through all documents
#' using the learned mu, sigma, and beta parameters. It optimizes the
#' document-level variational parameters (eta/lambda) for each document
#' and computes fresh theta values.
#'
#' Note: This does NOT re-optimize the global parameters (mu, sigma, beta),
#' only the document-level parameters.
#'
#' @examples
#' \dontrun{
#' # Fit large corpus without computing theta
#' model <- stm_svi(docs, vocab, K=50, compute_final_theta=FALSE)
#'
#' # Later, compute theta when needed
#' model <- update_svi_theta(model, docs, cores=4)
#'
#' # Now can use theta-dependent functions
#' findThoughts(model, n=5)
#' plot(model, type="summary")
#' }
#'
#' @seealso \code{\link{stm_svi}}
#'
#' @export
update_svi_theta <- function(model, documents=NULL, cores=NULL, verbose=TRUE) {
  # Validate input
  if(!inherits(model, "STM")) {
    stop("model must be an STM object")
  }

  if(is.null(model$svi) || !model$svi) {
    warning("Model does not appear to be fitted with SVI. ",
            "This function is designed for SVI models.")
  }

  # Extract documents if not provided
  if(is.null(documents)) {
    # Try to get from model (may not be stored)
    if(!is.null(model$documents)) {
      documents <- model$documents
    } else {
      stop("documents argument required (not stored in model)")
    }
  }

  # Extract parameters from model
  mu <- model$mu$mu
  sigma <- model$sigma

  # Extract beta (handle both formats)
  if(!is.null(model$beta$logbeta)) {
    # Convert from log space
    beta <- lapply(model$beta$logbeta, exp)
  } else if(!is.null(model$beta$beta)) {
    beta <- model$beta$beta
  } else {
    stop("Cannot extract beta from model")
  }

  # Extract betaindex
  if(!is.null(model$settings$covariates$betaindex)) {
    betaindex <- model$settings$covariates$betaindex
  } else {
    # Default: all documents use aspect 1
    betaindex <- rep(1, length(documents))
  }

  # Use current eta as warm start (fallback to zeros if not available)
  lambda_init <- model$eta
  if(is.null(lambda_init)) {
    K <- model$settings$dim$K
    lambda_init <- matrix(0, nrow=length(documents), ncol=K - 1)
  }

  # Determine parallel chunks
  N <- length(documents)
  parallel_chunks <- determine_parallel_chunks(cores, N)

  if(verbose) {
    cat("\n")
    cat("==========================================\n")
    cat("Updating theta via final E-step\n")
    cat("==========================================\n")
    cat(sprintf("Documents: %d\n", length(documents)))
    cat(sprintf("Topics: %d\n", model$settings$dim$K))
    if(parallel_chunks > 1) {
      cat(sprintf("Using %d parallel chunks\n", parallel_chunks))
    }
    cat("\n")
  }

  # Perform final E-step (parallel)
  t0 <- proc.time()
  final_result <- perform_final_estep(
    documents = documents,
    betaindex = betaindex,
    mu = mu,
    sigma = sigma,
    beta = beta,
    lambda_init = lambda_init,
    parallel_chunks = parallel_chunks,
    verbose = verbose
  )
  elapsed <- (proc.time() - t0)[3]

  if(verbose) {
    cat("==========================================\n")
    cat(sprintf("Final E-step completed in %.1fs\n", elapsed))
    cat(sprintf("Mean ELBO: %.2f\n", mean(final_result$bound)))
    cat("==========================================\n\n")
  }

  # Update model with new theta and eta
  model$theta <- final_result$theta
  model$eta <- final_result$lambda

  # Add metadata about update
  if(is.null(model$convergence$final_bound)) {
    model$convergence$final_bound <- sum(final_result$bound)
  }
  model$convergence$theta_updated <- TRUE
  model$convergence$theta_update_time <- elapsed

  return(model)
}
