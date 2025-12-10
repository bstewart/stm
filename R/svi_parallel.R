# Parallel Processing for SVI Models - Final E-Step Only
#
# This file provides parallel processing ONLY for final theta computation
# via update_svi_theta(). Mini-batch processing in stm_svi() remains serial.

#' Determine Optimal Parallel Chunks for Final E-Step
#'
#' Auto-detects available cores and determines optimal number of parallel
#' chunks based on workload size. Only used for final theta computation.
#'
#' @param cores User-specified cores (NULL for auto-detect, 1 for serial, N for N cores)
#' @param n_docs Number of documents to process
#' @param min_docs_per_core Minimum documents per core (default: 400)
#'
#' @return Integer number of parallel chunks to use
#'
#' @keywords internal
determine_parallel_chunks <- function(cores = NULL, n_docs = NULL,
                                      min_docs_per_core = 400) {
  if(is.null(cores)) {
    # Auto-detect with smart caps
    available_cores <- parallel::detectCores()
    if(is.na(available_cores)) available_cores <- 1

    # Use all cores minus 1, but cap at 8 to avoid excessive overhead
    suggested_cores <- min(max(1, available_cores - 1), 8)

    # Don't create more chunks than makes sense for workload
    if(!is.null(n_docs)) {
      max_useful_chunks <- max(1, floor(n_docs / min_docs_per_core))
      return(min(suggested_cores, max_useful_chunks))
    }

    return(suggested_cores)

  } else if(cores == 1) {
    return(1)  # Explicit serial

  } else if(cores > 1) {
    return(as.integer(cores))

  } else {
    stop("cores must be NULL (auto-detect), 1 (serial), or > 1 (parallel)")
  }
}


#' Combine Chunked Sufficient Statistics
#'
#' Combines results from parallel E-step chunks into single sufficient
#' statistics structure. Used after parallel::mclapply() returns list of
#' estep results.
#'
#' @param chunk_results List of estep() results from parallel chunks
#'
#' @return Combined sufficient statistics list
#'
#' @keywords internal
combine_chunked_suffstats <- function(chunk_results) {
  # Combine sigma sufficient statistics (matrix addition)
  combined_sigma <- Reduce("+", lapply(chunk_results, function(x) x$sigma))

  # Combine beta sufficient statistics (list of matrices)
  A <- length(chunk_results[[1]]$beta)
  combined_beta <- vector("list", length = A)
  for(a in 1:A) {
    combined_beta[[a]] <- Reduce("+", lapply(chunk_results,
                                              function(x) x$beta[[a]]))
  }

  # Combine bound (concatenate vectors)
  combined_bound <- unlist(lapply(chunk_results, function(x) x$bound))

  # Combine lambda (row-bind matrices)
  combined_lambda <- do.call(rbind, lapply(chunk_results, function(x) x$lambda))

  return(list(
    sigma = combined_sigma,
    beta = combined_beta,
    bound = combined_bound,
    lambda = combined_lambda
  ))
}


#' Parallel E-Step for Final Theta Computation Only
#'
#' Uses socket-based parallelization (parLapply) to avoid Rcpp deadlock issues.
#' Only beneficial for large N (>= 2000 docs) due to ~2s socket overhead.
#'
#' @param documents Document list
#' @param doc_indices Indices of documents to process
#' @param betaindex Beta index for documents
#' @param beta Beta parameters
#' @param lambda_old Old lambda values
#' @param mu Mu parameters (can be (K-1)×1 or (K-1)×N for document-specific priors)
#' @param sigma Sigma parameters
#' @param parallel_chunks Number of parallel chunks
#' @param verbose Print progress
#' @param update_mu Whether mu is document-specific (TRUE) or global (FALSE)
#'
#' @return Sufficient statistics list
#'
#' @keywords internal
estep_parallel <- function(documents, doc_indices, betaindex, beta,
                           lambda_old, mu, sigma, parallel_chunks = 1,
                           verbose = FALSE, update_mu = FALSE) {

  n_docs <- length(doc_indices)

  # Check if parallelization makes sense
  # Socket clusters have high overhead (~2s per call), so need larger batches
  min_docs_per_core <- 400  # Increased from 100 based on testing
  use_parallel <- (parallel_chunks > 1 &&
                   n_docs >= parallel_chunks * min_docs_per_core)

  if(!use_parallel) {
    # Fall back to sequential estep
    # Extract appropriate mu subset if document-specific
    if(update_mu) {
      mu_subset <- mu[, doc_indices, drop=FALSE]
    } else {
      mu_subset <- mu
    }

    estep_fn <- if(exists("estep", mode="function")) estep else stm:::estep
    return(estep_fn(
      documents = documents[doc_indices],
      beta.index = betaindex[doc_indices],
      update.mu = update_mu,
      beta = beta,
      lambda.old = lambda_old[doc_indices, , drop=FALSE],
      mu = mu_subset,
      sigma = sigma,
      verbose = verbose
    ))
  }

  # Partition documents into chunks
  chunk_size <- ceiling(n_docs / parallel_chunks)
  chunk_positions <- split(seq_len(n_docs),
                           ceiling(seq_len(n_docs) / chunk_size))

  # Create socket cluster (works with Rcpp, cross-platform)
  cl <- parallel::makeCluster(parallel_chunks, type="PSOCK")
  on.exit(parallel::stopCluster(cl), add=TRUE)

  # Export necessary objects to cluster nodes
  parallel::clusterExport(cl, c("documents", "doc_indices", "betaindex",
                                "beta", "lambda_old", "mu", "sigma",
                                "update_mu"),
                         envir=environment())

  # Load stm package on each node
  parallel::clusterEvalQ(cl, library(stm))

  # Process chunks in parallel using parLapply
  chunk_results <- parallel::parLapply(cl, chunk_positions, function(chunk_pos) {
    # Use stm:::estep (always available since we loaded stm on nodes)
    chunk_idx <- doc_indices[chunk_pos]

    # Extract appropriate mu for this chunk
    if(update_mu) {
      mu_chunk <- mu[, chunk_idx, drop=FALSE]
    } else {
      mu_chunk <- mu
    }

    stm:::estep(
      documents = documents[chunk_idx],
      beta.index = betaindex[chunk_idx],
      update.mu = update_mu,
      beta = beta,
      lambda.old = lambda_old[chunk_idx, , drop=FALSE],
      mu = mu_chunk,
      sigma = sigma,
      verbose = FALSE  # Suppress verbose in parallel
    )
  })

  # Combine results from chunks
  combined <- combine_chunked_suffstats(chunk_results)

  return(combined)
}
