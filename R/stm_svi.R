#' Compute Adaptive SVI Hyperparameters
#'
#' Internal function to automatically select reasonable SVI hyperparameters
#' based on corpus characteristics (N, K, V) and covariate presence.
#'
#' @param N Number of documents
#' @param K Number of topics
#' @param V Vocabulary size
#' @param has_prevalence Whether prevalence covariates are used
#' @return List with adaptive hyperparameter values
#'
#' @keywords internal
compute_svi_defaults <- function(N, K, V, has_prevalence) {
  # 1. Batch size (scales with topic complexity)
  base_batch <- max(64, min(512, 4 * K))
  if(N < 5000) {
    batch_size <- min(base_batch, max(32, floor(N * 0.2)))
  } else if(N < 50000) {
    batch_size <- base_batch
  } else {
    batch_size <- min(base_batch * 2, 1024)
  }

  # 2. Learning rate (adjusts for gradient scaling)
  scale_factor <- N / batch_size
  if(scale_factor < 50) {
    lr <- min(0.02, 0.01 * 1.5)
  } else if(scale_factor > 500) {
    lr <- max(0.005, 0.01 * 0.7)
  } else {
    lr <- 0.01
  }
  if(K > 75) lr <- lr * 0.8

  # 3. Max epochs (targets 500-1000 total iterations)
  iters_per_epoch <- ceiling(N / batch_size)
  if(N < 2000) {
    target_iters <- 1000
  } else if(N < 10000) {
    target_iters <- 800
  } else if(N < 50000) {
    target_iters <- 600
  } else {
    target_iters <- 500
  }
  max_epochs <- max(5, min(100, ceiling(target_iters / iters_per_epoch)))

  # 4. Convergence window (larger for smaller batches)
  convergence_window <- sqrt(128 / batch_size) * 10
  if(N > 50000) convergence_window <- convergence_window * 0.8
  convergence_window <- max(5, min(20, round(convergence_window)))

  # 5. Patience (higher for noisier small batches)
  patience <- 20 * sqrt(128 / batch_size)
  if(N < 2000) patience <- patience * 0.7
  patience <- max(ceiling(2 * convergence_window), round(patience))
  patience <- max(10, min(50, patience))

  # 6. Gamma update frequency (for prevalence covariates)
  if(has_prevalence) {
    if(N < 5000) {
      gamma_update_every <- max(5, ceiling(iters_per_epoch * 0.5))
    } else if(N < 50000) {
      gamma_update_every <- max(10, ceiling(iters_per_epoch))
    } else {
      gamma_update_every <- ceiling(iters_per_epoch * 2)
    }
    gamma_update_every <- max(5, min(100, gamma_update_every))
  } else {
    gamma_update_every <- NULL
  }

  return(list(
    batch_size = batch_size,
    lr = lr,
    max_epochs = max_epochs,
    patience = patience,
    convergence_window = convergence_window,
    gamma_update_every = gamma_update_every,
    iters_per_epoch = iters_per_epoch,
    target_total_iters = max_epochs * iters_per_epoch,
    gradient_scale = scale_factor
  ))
}

#' Stochastic Variational Inference for Structural Topic Model
#'
#' Estimates the Structural Topic Model using stochastic variational inference
#' with Adam optimizer. This enables efficient estimation on large corpora by
#' processing mini-batches of documents rather than the full dataset.
#'
#' This function implements stochastic variational EM with mini-batch updates
#' and adaptive learning rates via the Adam optimizer. It is particularly
#' suited for large document collections (N > 10,000) where standard EM becomes
#' computationally expensive.
#'
#' @param documents The document term matrix in STM format (list of matrices)
#' @param vocab Character vector of vocabulary terms
#' @param K Number of topics (must be >= 2)
#' @param prevalence Prevalence covariate formula (NOT YET SUPPORTED in prototype)
#' @param content Content covariate formula (NOT YET SUPPORTED in prototype)
#' @param data Data frame with covariates
#' @param init.type Initialization method: "Spectral" (default), "LDA", or "Random"
#' @param seed Random seed for reproducibility
#' @param batch_size Mini-batch size for stochastic updates. If NULL (default),
#'   automatically determined based on K and N (typically scales as 4×K, adjusted
#'   for corpus size). Reasonable range: 32-1024.
#' @param max_epochs Maximum number of passes through the data. If NULL (default),
#'   automatically targets 500-1000 total iterations based on N. Small corpora
#'   need more epochs; large corpora need fewer.
#' @param max_iters Alternative to max_epochs: maximum iterations (NULL by default)
#' @param lr Learning rate for Adam optimizer. If NULL (default), automatically
#'   determined based on N/batch_size gradient scaling and K (typically 0.005-0.02).
#'   Larger values converge faster but may be unstable.
#' @param adam_beta1 Adam first moment decay rate (default: 0.9)
#' @param adam_beta2 Adam second moment decay rate (default: 0.999)
#' @param adam_epsilon Adam numerical stability constant (default: 1e-8)
#' @param patience Number of iterations without improvement before early stopping.
#'   If NULL (default), automatically scales with batch_size (smaller batches need
#'   more patience due to noisier convergence). Typical range: 10-50.
#' @param convergence_window Number of recent iterations to average for convergence
#'   detection. If NULL (default), adapts to batch_size (larger window for smaller
#'   batches provides more smoothing). Typical range: 5-20.
#' @param eval_every Compute full ELBO every N iterations (NULL = never, expensive)
#' @param compute_final_theta Perform final E-step on all documents after convergence
#'   to compute accurate theta. If NULL (default), automatically set to TRUE for
#'   N < 10,000 and FALSE for N >= 10,000.
#'
#'   **Important**: If FALSE, the returned model will have \code{theta=NULL} and
#'   \code{eta=NULL}. Many downstream functions (findThoughts, plot.STM, topicCorr,
#'   etc.) require theta and will raise an error. To compute theta later, use:
#'   \code{model <- update_svi_theta(model, documents)}
#'
#'   Set to TRUE to always perform final E-step (slower but provides complete model),
#'   or FALSE to skip (faster, saves memory, but theta must be computed separately).
#' @param gamma_update_every How often to update gamma coefficients when using prevalence
#'   covariates (in iterations). If NULL (default), automatically adapts to corpus size:
#'   more frequent for small N (cheap), less frequent for large N (expensive). Gamma is
#'   updated by regressing lambda on covariates. Larger values = less frequent updates
#'   (faster but potentially less accurate).
#' @param verbose Print progress information
#' @param reportevery Report progress every N iterations
#' @param LDAbeta Use LDA-style beta (default: TRUE, required for prototype)
#' @param gamma.prior Prior for prevalence regression: "Pooled" (default, Bayesian
#'   regression with half-Cauchy priors) or "L1" (L1-regularized via glmnet).
#'   Only used when prevalence covariates are specified. Note: L1 mode is
#'   experimental; Pooled mode is recommended for most applications.
#' @param sigma.prior Strength of regularization toward diagonal covariance (0-1)
#' @param control Additional control parameters (list)
#' @param model Optional pre-fit STM model to resume from
#'
#' @return An STM object with additional $svi=TRUE flag
#'
#' @details
#' The stochastic variational inference algorithm processes documents in
#' mini-batches, computing gradients from sufficient statistics and updating
#' parameters via the Adam optimizer. The critical scaling factor N/batch_size
#' ensures unbiased gradient estimates.
#'
#' Key differences from standard \code{\link{stm}}:
#' \itemize{
#'   \item Uses mini-batch updates instead of full-batch EM
#'   \item Adaptive learning rates via Adam (no manual tuning needed)
#'   \item Convergence based on windowed ELBO and patience-based early stopping
#'   \item Much faster for large N (10k+ documents)
#'   \item May have slightly higher variance in final estimates
#' }
#'
#' @section Automatic Hyperparameter Selection:
#' When hyperparameters are left as NULL (default), they are automatically
#' selected based on corpus characteristics:
#'
#' \itemize{
#'   \item \strong{batch_size}: Scales with number of topics (4×K base), adjusted
#'     for corpus size. Small corpora (N<5K) use ≤20\% of documents to avoid waste.
#'     Large corpora (N>50K) can use larger batches (up to 2×base) for stability.
#'   \item \strong{lr}: Adjusted for gradient scaling (N/batch_size) and topic
#'     complexity. Automatically reduced for large corpora or when K>75 to prevent
#'     overshooting.
#'   \item \strong{max_epochs}: Targets 500-1000 total iterations depending on N.
#'     Small corpora need more passes through data; large corpora converge in fewer
#'     epochs.
#'   \item \strong{patience}: Scales with batch_size using sqrt relationship.
#'     Smaller batches produce noisier ELBO estimates and need higher patience
#'     to avoid premature stopping.
#'   \item \strong{convergence_window}: Smoothing window for ELBO, larger for
#'     smaller batches to reduce noise in convergence detection.
#'   \item \strong{gamma_update_every}: Adapts to corpus size. More frequent for
#'     small N (cheap gamma updates), less frequent for large N (expensive).
#' }
#'
#' These defaults work well for most applications. Manual tuning may improve
#' performance for specific use cases. Setting any parameter explicitly overrides
#' the automatic selection.
#'
#' @examples
#' \dontrun{
#' # Prepare corpus
#' data(poliblog5k)
#' out <- prepDocuments(poliblog5k.docs, poliblog5k.voc, poliblog5k.meta)
#'
#' # Fit with SVI
#' model_svi <- stm_svi(
#'   documents = out$documents,
#'   vocab = out$vocab,
#'   K = 20
#' )
#'
#' # Use like regular STM
#' labelTopics(model_svi)
#' plot(model_svi)
#' }
#'
#' @seealso \code{\link{stm}} for standard EM estimation
#'
#' @references
#' Hoffman et al. (2013) "Stochastic Variational Inference" JMLR
#' Kingma & Ba (2014) "Adam: A Method for Stochastic Optimization" ICLR
#'
#' @export
stm_svi <- function(documents, vocab, K,
                    prevalence=NULL, content=NULL, data=NULL,
                    init.type=c("Spectral", "LDA", "Random"),
                    seed=NULL,
                    batch_size=NULL,
                    max_epochs=NULL,
                    max_iters=NULL,
                    lr=NULL,
                    adam_beta1=0.9,
                    adam_beta2=0.999,
                    adam_epsilon=1e-8,
                    patience=NULL,
                    convergence_window=NULL,
                    eval_every=NULL,
                    compute_final_theta=NULL,
                    gamma_update_every=NULL,
                    verbose=TRUE,
                    reportevery=5,
                    LDAbeta=TRUE,
                    gamma.prior=c("Pooled", "L1"),
                    sigma.prior=0,
                    control=list(),
                    model=NULL) {

  # Match arguments and save call
  init.type <- match.arg(init.type)
  gamma.prior <- match.arg(gamma.prior)
  Call <- match.call()

  # Convert corpus to internal STM format
  args <- asSTMCorpus(documents, vocab, data)
  documents <- args$documents
  vocab <- args$vocab
  data <- args$data

  # --- Validate documents ---
  if(missing(documents)) stop("Must include documents")
  if(!is.list(documents)) stop("documents must be a list, see documentation.")
  if(!all(unlist(lapply(documents, is.matrix)))) {
    stop("Each list element in documents must be a matrix. See documentation.")
  }
  if(any(unlist(lapply(documents, function(x) anyDuplicated(x[1,]))))) {
    stop("Duplicate term indices within a document. See documentation.")
  }
  N <- length(documents)

  # --- Extract and check word indices ---
  wcountvec <- unlist(lapply(documents, function(x) rep(x[1,], times=x[2,])),
                     use.names=FALSE)
  wcounts <- list(Group.1=sort(unique(wcountvec)))
  V <- length(wcounts$Group.1)

  if(!posint(wcounts$Group.1)) {
    stop("Word indices are not positive integers")
  }
  if(!isTRUE(all.equal(wcounts$Group.1, 1:V))) {
    stop("Word indices must be sequential integers starting with 1.")
  }
  wcounts$x <- tabulate(wcountvec)
  rm(wcountvec)

  # --- Check vocabulary ---
  if(length(vocab) != V) {
    stop("Vocab length does not match observed word indices")
  }

  # --- Check number of topics ---
  if(missing(K)) stop("K, the number of topics, is required.")
  if(K != 0) {
    if(!(posint(K) && length(K)==1 && K>1)) {
      stop("K must be a positive integer greater than 1.")
    }
    if(K == 2) {
      warning("K=2 is equivalent to a unidimensional scaling model.")
    }
  } else {
    stop("K=0 (automatic topic selection) not yet supported in stm_svi")
  }

  # --- Compute adaptive defaults if parameters not specified ---
  defaults <- compute_svi_defaults(
    N = N,
    K = K,
    V = V,
    has_prevalence = !is.null(prevalence)
  )

  # Apply defaults only if user didn't specify
  if(is.null(batch_size)) {
    batch_size <- defaults$batch_size
    if(verbose) {
      cat(sprintf("Using adaptive batch_size=%d (for K=%d, N=%d)\n",
                  batch_size, K, N))
    }
  }

  if(is.null(lr)) {
    lr <- defaults$lr
    if(verbose) {
      cat(sprintf("Using adaptive lr=%.4f (gradient scale=%.1f)\n",
                  lr, defaults$gradient_scale))
    }
  }

  if(is.null(max_epochs)) {
    max_epochs <- defaults$max_epochs
    if(verbose) {
      cat(sprintf("Using adaptive max_epochs=%d (~%d total iters)\n",
                  max_epochs, defaults$target_total_iters))
    }
  }

  if(is.null(patience)) {
    patience <- defaults$patience
  }

  if(is.null(convergence_window)) {
    convergence_window <- defaults$convergence_window
  }

  if(is.null(gamma_update_every) && !is.null(prevalence)) {
    gamma_update_every <- defaults$gamma_update_every
  }

  # Validation: batch_size can't exceed N
  if(batch_size > N) {
    batch_size <- N
    if(verbose) {
      cat(sprintf("Adjusting batch_size to %d (corpus size)\n", N))
    }
  }

  # --- Validate SVI-specific parameters ---
  if(!is.numeric(batch_size) || batch_size < 1) {
    stop("batch_size must be a positive integer")
  }
  if(batch_size > N) {
    stop(sprintf("batch_size (%d) must be <= number of documents (%d)", batch_size, N))
  }
  if(batch_size < 10) {
    warning("Very small batch_size may lead to unstable optimization")
  }

  if(!is.numeric(lr) || lr <= 0) {
    stop("Learning rate (lr) must be positive")
  }
  if(lr > 0.1) {
    warning("Large learning rate may cause divergence. Consider lr < 0.1")
  }

  if(!is.null(max_epochs) && (!is.numeric(max_epochs) || max_epochs < 1)) {
    stop("max_epochs must be a positive number")
  }

  if(!is.numeric(patience) || patience < 1) {
    stop("patience must be a positive integer")
  }

  # --- Process prevalence covariates ---
  if(!is.null(prevalence)) {
    # Use makeTopMatrix to process formula
    covariates <- list()
    covariates$X <- makeTopMatrix(prevalence, data)

    if(ncol(covariates$X) == 1) {
      # Only intercept - equivalent to no covariates
      if(verbose) {
        cat("Note: Prevalence formula contains only intercept, treating as CTM mode\n")
      }
      prevalence <- NULL
      covariates <- NULL
    }
  } else {
    covariates <- NULL
  }

  # --- Check content covariates (not yet supported) ---
  if(!is.null(content)) {
    stop("Content covariates not yet supported in stm_svi prototype. ",
         "Use standard stm() function or set content=NULL.")
  }

  # For prototype: force simple CTM mode (no content covariates)
  betaindex <- rep(1, N)
  A <- 1

  # --- Check other parameters ---
  if(!is.logical(verbose)) stop("verbose must be logical")
  if(!is.logical(LDAbeta)) stop("LDAbeta must be logical")
  if(!LDAbeta) {
    stop("LDAbeta=FALSE not yet supported in stm_svi prototype")
  }
  if(sigma.prior < 0 || sigma.prior > 1) {
    stop("sigma.prior must be between 0 and 1")
  }

  # --- Process seed ---
  if(is.null(seed)) {
    seed <- floor(runif(1) * 1e7)
  }
  set.seed(seed)

  # --- Process compute_final_theta parameter ---
  if(is.null(compute_final_theta)) {
    # Default: TRUE for small corpora, FALSE for large
    compute_final_theta <- (N < 10000)
    if(verbose) {
      cat(sprintf("Auto-setting compute_final_theta=%s (N=%d)\n",
                  compute_final_theta, N))
    }
  }
  if(!is.logical(compute_final_theta)) {
    stop("compute_final_theta must be TRUE, FALSE, or NULL (auto)")
  }

  # --- Construct settings object ---
  settings <- list(
    dim = list(K=K, A=A, V=V, N=N, wcounts=wcounts),
    verbose = verbose,
    topicreportevery = reportevery,

    # Standard convergence settings (mostly ignored in SVI)
    convergence = list(
      max.em.its = NULL,
      em.converge.thresh = NULL,
      allow.neg.change = TRUE
    ),

    # Covariate settings
    covariates = list(
      X = if(!is.null(covariates)) covariates$X else NULL,
      betaindex = betaindex,
      yvarlevels = NULL,
      formula = prevalence  # Store the formula object for fitNewDocuments
    ),

    # Gamma settings (prevalence prior)
    gamma = list(
      mode = if(is.null(covariates)) "CTM" else gamma.prior,
      prior = NULL,
      enet = 1,      # For L1 mode: 1=lasso, 0=ridge, 0<x<1=elastic net
      ic.k = 2,      # For L1 mode: IC selection (2=BIC-like)
      maxits = 1000,
      update_every = gamma_update_every
    ),

    # Sigma settings (covariance prior)
    sigma = list(prior=sigma.prior),

    # Kappa settings (content prior - using LDA mode)
    kappa = list(
      LDAbeta = LDAbeta,
      interactions = FALSE,
      fixedintercept = TRUE,
      mstep = list(tol=.001, maxit=3),
      contrast = FALSE
    ),

    # Tau settings (content covariate prior)
    tau = list(
      mode = "L1",
      tol = 1e-5,
      enet = 1,
      nlambda = 250,
      lambda.min.ratio = .001,
      ic.k = 2,
      maxit = 1e4
    ),

    # Initialization settings (must match stm.init expectations)
    init = list(
      mode = init.type,
      nits = 50,
      burnin = 25,
      alpha = 50/K,
      eta = 0.01,
      s = .05,
      p = 3000,
      d.group.size = 2000,
      recoverEG = TRUE,
      tSNE_init.dims = 50,
      tSNE_perplexity = 30
    ),

    # SVI-specific settings
    svi = list(
      batch_size = batch_size,
      max_epochs = max_epochs,
      max_iters = max_iters,
      lr = lr,
      adam_beta1 = adam_beta1,
      adam_beta2 = adam_beta2,
      adam_epsilon = adam_epsilon,
      patience = patience,
      convergence_window = convergence_window,
      eval_every = eval_every,
      compute_final_theta = compute_final_theta
    ),

    seed = seed,
    ngroups = 1,  # For parallelization (not used in SVI)
    call = Call,
    prevalence = covariates  # Store design matrix for prevalence covariates
  )

  # Add spectral maxV if needed
  if(init.type == "Spectral" && V > 10000) {
    settings$init$maxV <- 10000
  }

  # --- Print startup message ---
  if(verbose) {
    cat("\n")
    cat("===================================================\n")
    cat("Stochastic Variational Inference for STM\n")
    cat("===================================================\n")
    cat(sprintf("Documents: %d, Vocabulary: %d, Topics: %d\n", N, V, K))
    cat(sprintf("Batch size: %d, Learning rate: %.4f\n", batch_size, lr))
    cat(sprintf("Max epochs: %d, Patience: %d\n", max_epochs, patience))
    cat(sprintf("Initialization: %s\n", init.type))
    if(!is.null(covariates)) {
      cat(sprintf("Prevalence: %d covariates (mode=%s, update every %d iters)\n",
                  ncol(covariates$X)-1, gamma.prior, gamma_update_every))
    } else {
      cat("Prevalence: None (CTM mode)\n")
    }
    cat("===================================================\n\n")
  }

  # --- Run SVI ---
  result <- stm.svi.control(documents, vocab, settings, model)

  return(result)
}
