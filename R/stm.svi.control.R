# Stochastic Variational Inference Control Loop for STM
#
# This is the main control loop for SVI, analogous to stm.control.R
# but implementing mini-batch stochastic updates with Adam optimizer.

#' SVI Control Loop
#'
#' Implements stochastic variational inference for STM using mini-batch
#' updates and Adam optimization. This is the core function that orchestrates
#' the SVI algorithm.
#'
#' @param documents Document list in STM format
#' @param vocab Vocabulary vector
#' @param settings Settings list from stm_svi()
#' @param model Optional: pre-fit model to resume from
#'
#' @return STM model object
#'
#' @keywords internal
stm.svi.control <- function(documents, vocab, settings, model=NULL) {
  # Extract dimensions and parameters
  N <- settings$dim$N
  K <- settings$dim$K
  A <- settings$dim$A
  V <- settings$dim$V
  betaindex <- settings$covariates$betaindex

  # Extract SVI-specific settings
  batch_size <- settings$svi$batch_size
  lr <- settings$svi$lr
  adam_beta1 <- settings$svi$adam_beta1
  adam_beta2 <- settings$svi$adam_beta2
  adam_epsilon <- settings$svi$adam_epsilon
  verbose <- settings$verbose

  # --- Step 1: Initialize model parameters ---
  if(verbose) cat("Initializing model...\n")

  if(is.null(model)) {
    # Use standard STM initialization
    # Use ::: to access internal function when sourcing files
    init_result <- if(exists("stm.init", mode="function")) {
      stm.init(documents, settings)
    } else {
      stm:::stm.init(documents, settings)
    }
    mu <- init_result$mu
    sigma <- init_result$sigma
    beta <- init_result$beta
    lambda <- init_result$lambda
    gamma <- init_result$gamma  # Will be NULL from standard init

    # If prevalence covariates are present, initialize gamma here
    if(!is.null(settings$prevalence)) {
      if(verbose) cat("  Initializing gamma for prevalence covariates...\n")
      X <- settings$prevalence$X
      # Initialize gamma with small random values
      gamma <- matrix(rnorm(ncol(X) * (K-1), mean=0, sd=0.1),
                      nrow=ncol(X), ncol=K-1)
      # Compute document-specific mu from gamma
      # X might be sparse, so ensure result is dense matrix
      mu_temp <- as.matrix(X %*% gamma)  # N × (K-1)
      mu <- t(mu_temp)  # (K-1) × N matrix

      if(verbose) {
        cat(sprintf("  Gamma: %d x %d, Mu: %d x %d\n",
                    nrow(gamma), ncol(gamma), nrow(mu), ncol(mu)))
      }
    }
  } else {
    # Resume from previous model
    mu <- if(!is.null(model$mu$mu)) model$mu$mu else model$mu
    sigma <- model$sigma
    if(!is.null(model$beta$logbeta)) {
      beta <- lapply(model$beta$logbeta, exp)
    } else if(!is.null(model$beta$beta)) {
      beta <- model$beta$beta
    } else {
      stop("Cannot extract beta from model for SVI resume")
    }
    lambda <- model$eta  # eta contains lambda
    gamma <- if(!is.null(model$mu$gamma)) model$mu$gamma else NULL
  }

  # --- Holdout setup (optional) ---
  holdout_idx <- integer(0)
  holdout_tokens <- NULL
  if(!is.null(settings$svi$eval_holdout_size)) {
    holdout_size <- as.integer(settings$svi$eval_holdout_size)
    holdout_seed <- settings$svi$eval_holdout_seed
    if(!is.null(holdout_seed)) {
      has_seed <- exists(".Random.seed", envir=.GlobalEnv)
      if(has_seed) old_seed <- .Random.seed
      set.seed(holdout_seed)
      holdout_idx <- sample(1:N, holdout_size, replace=FALSE)
      if(has_seed) {
        .Random.seed <<- old_seed
      } else {
        rm(.Random.seed, envir=.GlobalEnv)
      }
    } else {
      holdout_idx <- sample(1:N, holdout_size, replace=FALSE)
    }
  }

  train_idx <- setdiff(1:N, holdout_idx)
  N_train <- length(train_idx)

  if(length(holdout_idx) > 0 && verbose) {
    cat(sprintf("Using %d holdout documents for ELBO monitoring\n", length(holdout_idx)))
  }
  if(length(holdout_idx) > 0) {
    holdout_tokens <- sum(vapply(documents[holdout_idx], function(x) sum(x[2, ]),
                                 numeric(1)))
  }
  train_tokens <- sum(vapply(documents[train_idx], function(x) sum(x[2, ]),
                             numeric(1)))
  settings$svi$train_tokens <- train_tokens

  # --- Step 2: Initialize Adam optimizer state ---
  if(verbose) cat("Initializing Adam optimizer...\n")
  has_prevalence <- !is.null(settings$prevalence)
  adam_state <- initialize_adam_state(mu, sigma, beta, has_prevalence=has_prevalence)

  # --- Step 3: Initialize convergence tracking ---
  convergence <- initialize_svi_convergence(settings)
  t0 <- proc.time()

  # Make beta and vocab accessible for progress reporting
  settings$beta_current <- NULL  # Will be updated during iterations
  settings$vocab <- vocab  # Store vocab in settings for reporting

  # --- Step 4: Main SVI loop ---
  if(batch_size > N_train) {
    batch_size <- N_train
  }

  if(verbose) {
    cat(sprintf("Beginning Stochastic Variational Inference (batch_size=%d, lr=%.4f)\n",
               batch_size, lr))
    cat(sprintf("  Documents: %d, Vocabulary: %d, Topics: %d\n", N_train, V, K))
    cat(sprintf("  Max epochs: %d (~%d iterations), Patience: %d\n",
               settings$svi$max_epochs, ceiling(N_train / batch_size) * settings$svi$max_epochs,
               settings$svi$patience))
  }

  iter <- 1
  while(!convergence$stopits) {
    # --- Sample mini-batch ---
    batch_idx <- sample(train_idx, min(batch_size, N_train), replace=FALSE)
    batch_docs <- documents[batch_idx]

    # --- Extract document-specific mu if covariates present ---
    if(!is.null(settings$prevalence)) {
      mu_batch <- mu[, batch_idx, drop=FALSE]  # (K-1) × batch_size
      update_mu_flag <- TRUE
    } else {
      mu_batch <- mu  # Scalar or (K-1) × 1 matrix
      update_mu_flag <- FALSE
    }

    # --- E-step on mini-batch (serial) ---
    suffstats <- tryCatch({
      estep_fn <- if(exists("estep", mode="function")) estep else stm:::estep
      estep_fn(
        documents = batch_docs,
        beta.index = betaindex[batch_idx],
        update.mu = update_mu_flag,
        beta = beta,
        lambda.old = lambda[batch_idx, , drop=FALSE],
        mu = mu_batch,
        sigma = sigma,
        verbose = FALSE
      )
    }, error = function(e) {
      cat(sprintf("\nError in E-step at iteration %d:\n", iter))
      cat(conditionMessage(e), "\n")
      stop("E-step failed. Check parameters for numerical issues.")
    })

    # --- Convert sufficient statistics to gradients ---
    grad_mu_input <- if(!is.null(settings$prevalence)) mu_batch else mu
    gradients <- compute_svi_gradients(
      suffstats, N_train, batch_size, grad_mu_input, sigma, beta, settings
    )

    # Check for numerical issues
    if(!check_gradients(gradients, iter, stop_on_error=FALSE)) {
      # If gradients are bad, try to recover
      cat("Warning: Bad gradients detected, clipping...\n")
      gradients <- clip_gradients(gradients, clip_value=5.0)
    }

    # --- Adam optimizer step ---
    # Apply learning rate decay to prevent late-stage instability
    # Without decay, constant LR causes parameters to drift as noise dominates
    # Use power-law decay following Robbins-Monro / Hoffman et al. (2013):
    #   lr_t = lr * (tau + 1) / (tau + t)^kappa
    # This satisfies convergence conditions: sum(lr_t) = inf, sum(lr_t^2) < inf
    lr_kappa <- if(!is.null(settings$svi$lr_kappa)) settings$svi$lr_kappa else 0.6
    lr_tau <- if(!is.null(settings$svi$lr_tau)) settings$svi$lr_tau else 10
    lr_decayed <- lr * (lr_tau + 1) / (lr_tau + iter)^lr_kappa

    adam_state <- adam_update_step(
      adam_state, gradients, lr_decayed, iter,
      adam_beta1, adam_beta2, adam_epsilon
    )

    # --- Apply parameter updates ---
    # Update mu (skip if NULL - prevalence covariates case)
    # Use gradient ASCENT (addition) for ELBO maximization
    if(!is.null(adam_state$update_mu)) {
      mu <- mu + adam_state$update_mu
    }

    # Update sigma
    # Use gradient ASCENT (addition) for ELBO maximization
    sigma <- sigma + adam_state$update_sigma
    sigma <- ensure_sigma_pd(sigma)  # Project back to PD cone if needed

    # Update beta (beta is a list of matrices)
    beta <- update_beta_from_adam(beta, adam_state$update_beta)

    # --- Update lambda for resampled documents ---
    # Use the new lambda values from E-step as warm starts
    lambda[batch_idx, ] <- suffstats$lambda

    # --- Convergence check ---
    convergence <- svi_convergence_check(
      sum(suffstats$bound), convergence, batch_size, N_train, settings
    )

    # --- Periodic gamma update (if covariates present) ---
    if(!is.null(settings$prevalence) &&
       iter %% settings$gamma$update_every == 0) {

      if(verbose) {
        cat(sprintf("  [Iteration %d: Updating gamma coefficients]\n", iter))
      }

      # Perform full E-step to get current lambda for all documents
      full_suffstats <- tryCatch({
        estep_fn(
          documents = documents,
          beta.index = betaindex,
          update.mu = TRUE,
          beta = beta,
          lambda.old = lambda,
          mu = mu,
          sigma = sigma,
          verbose = FALSE
        )
      }, error = function(e) {
        cat(sprintf("Warning: Full E-step failed during gamma update at iteration %d\n", iter))
        cat(conditionMessage(e), "\n")
        return(NULL)
      })

      if(!is.null(full_suffstats)) {
        # Update lambda
        lambda <- full_suffstats$lambda

        # Update gamma via regression
        opt.mu_fn <- if(exists("opt.mu", mode="function")) opt.mu else stm:::opt.mu

        mu_result <- tryCatch({
          opt.mu_fn(
            lambda = lambda,
            mode = settings$gamma$mode,
            covar = settings$prevalence$X,
            enet = settings$gamma$enet,
            ic.k = settings$gamma$ic.k,
            maxits = settings$gamma$maxits
          )
        }, error = function(e) {
          cat(sprintf("Warning: opt.mu failed at iteration %d\n", iter))
          cat(conditionMessage(e), "\n")
          return(NULL)
        })

        if(!is.null(mu_result)) {
          # Update mu and gamma
          mu <- mu_result$mu  # (K-1) × N matrix
          gamma <- mu_result$gamma  # (P+1) × (K-1) matrix

          if(verbose) {
            cat(sprintf("    Gamma updated (mode=%s)\n", settings$gamma$mode))
          }
        }
      }
    }

    # --- Report progress ---
    # Update beta reference for reporting
    settings$beta_current <- beta

    report_svi_progress(iter, convergence, settings)

    # --- Periodic holdout ELBO evaluation ---
    holdout_every <- settings$svi$eval_holdout_every
    if(length(holdout_idx) > 0 && !is.null(holdout_every) &&
       iter %% holdout_every == 0) {
      estep_fn <- if(exists("estep", mode="function")) estep else stm:::estep
      update_mu_holdout <- !is.null(settings$prevalence)
      mu_holdout <- if(update_mu_holdout) mu[, holdout_idx, drop=FALSE] else mu
      holdout_suffstats <- estep_fn(
        documents = documents[holdout_idx],
        beta.index = betaindex[holdout_idx],
        update.mu = update_mu_holdout,
        beta = beta,
        lambda.old = lambda[holdout_idx, , drop=FALSE],
        mu = mu_holdout,
        sigma = sigma,
        verbose = FALSE
      )

      lambda[holdout_idx, ] <- holdout_suffstats$lambda
      holdout_bound <- sum(holdout_suffstats$bound)
      holdout_per_token <- if(!is.null(holdout_tokens) && holdout_tokens > 0) {
        holdout_bound / holdout_tokens
      } else {
        NA_real_
      }
      convergence$holdout_bound_history <- c(convergence$holdout_bound_history,
                                             holdout_per_token)

      if(verbose) {
        if(is.na(holdout_per_token)) {
          cat(sprintf("  Holdout ELBO: %.2f\n", holdout_bound))
        } else {
          cat(sprintf("  Holdout ELBO (per token): %.4f\n", holdout_per_token))
        }
      }

      # Update holdout window and patience
      window_size <- if(!is.null(settings$svi$convergence_window)) {
        settings$svi$convergence_window
      } else {
        10
      }
      if(length(convergence$holdout_bound_history) >= window_size &&
         !is.na(holdout_per_token)) {
        recent <- tail(convergence$holdout_bound_history, window_size)
        window_mean <- mean(recent, na.rm=TRUE)
        if(!is.finite(window_mean)) {
          next
        }
        convergence$holdout_window <- c(convergence$holdout_window, window_mean)

        improvement_tol <- if(!is.null(settings$svi$improvement_tol)) {
          settings$svi$improvement_tol
        } else {
          1e-4
        }
        if(!is.finite(improvement_tol)) {
          improvement_tol <- 1e-4
        }
        if(!is.finite(convergence$holdout_best)) {
          convergence$holdout_best <- -Inf
        }
        improvement_threshold <- improvement_tol * max(1, abs(convergence$holdout_best))
        if(!is.finite(improvement_threshold)) {
          improvement_threshold <- 0
        }

        if(!is.finite(convergence$holdout_best)) {
          improved <- is.finite(window_mean)
        } else {
          improved <- isTRUE(
            is.finite(window_mean) && is.finite(improvement_threshold) &&
            window_mean > convergence$holdout_best + improvement_threshold
          )
        }
        if(improved) {
          convergence$holdout_best <- window_mean
          convergence$holdout_patience_counter <- 0
        } else {
          convergence$holdout_patience_counter <- convergence$holdout_patience_counter + 1
        }

        holdout_patience <- if(!is.null(settings$svi$holdout_patience)) {
          settings$svi$holdout_patience
        } else {
          settings$svi$patience
        }
        if(!is.finite(holdout_patience)) {
          holdout_patience <- settings$svi$patience
        }
        holdout_patience_iters <- holdout_patience
        if(!is.null(holdout_every)) {
          holdout_patience_iters <- holdout_patience / holdout_every
        }
        if(verbose) {
          cat(sprintf(
            "  Holdout window mean: %.4f (best: %.4f, tol: %.4g, patience: %.0f/%.0f)\n",
            window_mean,
            convergence$holdout_best,
            improvement_threshold,
            convergence$holdout_patience_counter,
            holdout_patience_iters
          ))
        }
        if(convergence$holdout_patience_counter >= holdout_patience_iters) {
          convergence$holdout_stop <- TRUE
          if(settings$svi$early_stop == "holdout") {
            convergence$converged <- TRUE
            convergence$stopits <- TRUE
            if(verbose) {
              cat(sprintf("\nConverged: holdout no improvement for %d evals (%.1f epochs)\n",
                         holdout_patience, convergence$epoch))
              cat(sprintf("Final holdout ELBO (per token): %.4f\n", window_mean))
            }
          }
        }
      }

      if(settings$svi$early_stop == "both" &&
         convergence$train_stop && convergence$holdout_stop) {
        convergence$converged <- TRUE
        convergence$stopits <- TRUE
        if(verbose) {
          cat(sprintf("\nConverged: train and holdout plateaued (%.1f epochs)\n",
                     convergence$epoch))
        }
      }
    }

    # --- Periodic full ELBO evaluation (optional, expensive) ---
    eval_every <- settings$svi$eval_every
    if(!is.null(eval_every) && iter %% eval_every == 0) {
      if(verbose) cat("  Computing full ELBO...")
      full_elbo <- compute_full_elbo(documents, mu, sigma, beta, lambda,
                                     betaindex, settings)
      if(verbose) cat(sprintf(" %.2f\n", full_elbo))
    }

    iter <- iter + 1

    # Safety check: prevent infinite loops
    if(iter > 100000) {
      warning("SVI exceeded 100,000 iterations without converging. Stopping.")
      convergence$stopits <- TRUE
    }
  }

  # --- Step 5: Compute final theta ---
  if(settings$svi$compute_final_theta) {
    if(verbose) cat("\nPerforming final E-step on all documents...\n")

    # Final gamma update if covariates present
    if(!is.null(settings$prevalence)) {
      if(verbose) cat("  Performing final gamma update...\n")

      opt.mu_fn <- if(exists("opt.mu", mode="function")) opt.mu else stm:::opt.mu

      mu_result <- tryCatch({
        opt.mu_fn(
          lambda = lambda,
          mode = settings$gamma$mode,
          covar = settings$prevalence$X,
          enet = settings$gamma$enet,
          ic.k = settings$gamma$ic.k,
          maxits = settings$gamma$maxits
        )
      }, error = function(e) {
        cat("Warning: Final gamma update failed, using last gamma values\n")
        cat(conditionMessage(e), "\n")
        NULL
      })

      if(!is.null(mu_result)) {
        mu <- mu_result$mu
        gamma <- mu_result$gamma
        if(verbose) cat("  Final gamma update complete.\n")
      }
    }

    # Perform complete E-step with final parameters (serial in stm_svi context)
    final_result <- perform_final_estep(
      documents = documents,
      betaindex = betaindex,
      mu = mu,
      sigma = sigma,
      beta = beta,
      lambda_init = lambda,
      parallel_chunks = 1,  # Force serial (parallelization only in update_svi_theta)
      verbose = verbose
    )

    # Use fresh lambda and theta
    lambda <- final_result$lambda
    theta <- final_result$theta

    if(verbose) {
      cat(sprintf("Final E-step complete. Mean ELBO: %.2f\n",
                  mean(final_result$bound)))
    }

  } else {
    # Don't compute theta - leave as NULL
    theta <- NULL
    lambda <- NULL

    if(verbose) {
      cat("\nSkipping final theta computation (compute_final_theta=FALSE)\n")
      cat("To compute theta later, use: update_svi_theta(model, documents)\n")
    }
  }

  # --- Step 6: Construct output STM object ---
  settings$svi$holdout_idx <- holdout_idx
  if(verbose) {
    if(convergence$converged) {
      cat("Model Converged \n")
    } else {
      cat("Model Terminated Before Convergence Reached \n")
    }
    cat(sprintf("Completed in %d iterations (%.1f epochs) \n",
               convergence$its, convergence$epoch))
  }

  # Convert beta to log space (STM standard format)
  # Create beta structure matching stm.control.R
  beta_list <- list()
  beta_list$logbeta <- vector("list", length(beta))
  for(i in 1:length(beta)) {
    beta_list$logbeta[[i]] <- log(beta[[i]])
  }

  # Calculate inverse sigma
  invsigma <- solve(sigma)

  # Record timing
  elapsed <- (proc.time() - t0)[3]

  # Construct output object matching STM format
  out <- list(
    mu = list(mu=mu, gamma=gamma),  # Wrap mu and gamma in list like stm.control.R
    sigma = sigma,
    beta = beta_list,  # Use beta list with $logbeta
    settings = settings,
    vocab = vocab,
    convergence = list(
      bound = convergence$bound_history,
      bound_window = convergence$bound_window,
      holdout_bound = convergence$holdout_bound_history,
      holdout_window = convergence$holdout_window,
      its = convergence$its,
      converged = convergence$converged,
      theta_computed = !is.null(theta)  # Track theta availability
    ),
    theta = theta,
    eta = lambda,
    invsigma = invsigma,
    time = elapsed,
    version = utils::packageVersion("stm"),
    svi = TRUE  # Flag indicating this was fit with SVI
  )

  class(out) <- "STM"
  return(out)
}
