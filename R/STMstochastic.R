# Stochastic variational EM for the structural topic model.
#
# Phase 1 is online EM (Cappe and Moulines 2009): the corpus is visited in
# shuffled minibatches, each minibatch's expected sufficient statistics are
# scaled to the corpus and blended into running estimates with a
# Robbins-Monro step size rho_t = (t + tau)^(-kappa), and the ordinary M-step
# is applied to the running estimates after every minibatch.  For the
# exponential-family M-steps used here this is the natural-gradient update of
# stochastic variational inference (Hoffman et al. 2013).  The phase ends when
# the smoothed per-token bound stops improving or after a fixed number of
# passes.
#
# Phase 2 hands the global parameters to the memoized incremental method
# (incremental.control), whose exact cached statistics converge to a
# stationary point of the batch objective.  Phase 1 does the fast approach
# from the initialization; phase 2 guarantees a solution of batch quality.

stochastic.defaults <- function(N, K, V) {
  # at least 20 minibatches per pass, never more than 256 documents each
  batch.size <- as.integer(max(64, min(256, ceiling(N/20))))
  list(batch.size=batch.size, kappa=0.7, tau=50, max.epochs=3L, tol=2e-3,
       finish=TRUE, checkpoints=4L, min.checkpoints=2L, patience=2L)
}


stochastic.control <- function(documents, vocab, settings) {
  globaltime <- proc.time()
  N <- settings$dim$N
  K <- settings$dim$K
  V <- settings$dim$V
  verbose <- settings$verbose
  betaindex <- settings$covariates$betaindex
  X <- settings$covariates$X
  control <- settings$stochastic

  if(verbose) cat(switch(EXPR=settings$init$mode,
    Spectral="Beginning Spectral Initialization \n",
    LDA="Beginning LDA Initialization \n",
    Random="Beginning Random Initialization \n",
    Custom="Beginning Custom Initialization \n"))
  initialized <- stm.init(documents, settings)
  mu <- list(mu=initialized$mu, gamma=NULL)
  sigma <- initialized$sigma
  beta <- initialized$beta
  lambda <- initialized$lambda
  rm(initialized)

  blocks <- incremental.blocks(N, control$batch.size, control$seed)
  nblocks <- length(blocks)
  record <- list(minibatch.bound=numeric(), minibatch.epoch=numeric(),
                 checkpoint.bound=numeric(), checkpoint.epoch=numeric(),
                 checkpoint.change=numeric(), epochs=0, updates=0L,
                 batch.size=control$batch.size, blocks=nblocks,
                 kappa=control$kappa, tau=control$tau,
                 stop.reason="skipped")

  if(nblocks < 8 || control$max.epochs == 0) {
    if(verbose) cat("Too few minibatches for a stochastic phase; using memoized inference directly.\n")
    record$stop.reason <- "too few minibatches"
    return(stochastic.finish(documents, vocab, settings, mu, sigma, beta,
                             lambda, record, globaltime))
  }

  has.covariates <- !is.null(X)
  ntokens <- sum(settings$dim$wcounts$x)
  block.tokens <- vapply(blocks, function(index) {
    sum(vapply(documents[index], function(document) sum(document[2,]),
               numeric(1)))
  }, numeric(1))

  # Running estimates of the complete-data sufficient statistics on the scale
  # of the whole corpus.  They start from the initialization: topic-word
  # counts implied by beta with the corpus mass split evenly across topics,
  # and the diffuse initial covariance.
  S.beta <- beta[[1]] * (ntokens/K)
  S.nu <- N * sigma
  S.lambda.sum <- numeric(K - 1)
  S.lambda.cross <- matrix(0, nrow=K - 1, ncol=K - 1)
  S.x.lambda <- if(has.covariates) matrix(0, nrow=ncol(X), ncol=K - 1) else NULL
  x.cross <- if(has.covariates) crossprod(X) else NULL
  # A small amount of mass spread evenly over the vocabulary keeps words that
  # have not yet been seen in a minibatch from getting probability zero.  It
  # is only used in this phase; the memoized phase uses exact statistics.
  smoothing <- .001 * ntokens/(K * V)

  checkpoint.every <- max(2L, as.integer(round(nblocks/control$checkpoints)))
  cache <- vector("list", nblocks)   # each block's statistics at its last visit
  bound <- rep(NA_real_, N)
  window <- numeric()
  window.gain <- numeric()
  block.bound <- rep(NA_real_, nblocks) # per-token bound of each block at its last visit
  ncheckpoints <- 0L
  flat <- 0L
  updates <- 0L
  stop.now <- FALSE
  stop.reason <- "maximum stochastic epochs reached"
  if(verbose) {
    cat(sprintf("Beginning stochastic inference with minibatches of %d documents (%d per pass).\n",
                control$batch.size, nblocks))
  }

  for(epoch in seq_len(control$max.epochs)) {
    block.order <- incremental.epoch.order(nblocks, control$seed, epoch)
    for(position in seq_along(block.order)) {
      b <- block.order[position]
      index <- blocks[[b]]
      # Before the first M-step the prior mean is the initialization's single
      # vector; afterwards it is the prevalence regression's prediction.
      update.mu <- !is.null(mu$gamma)
      mu.block <- if(update.mu) {
        t(as.matrix(X[index,,drop=FALSE] %*% mu$gamma))
      } else {
        mu$mu
      }
      ss <- estep(documents[index], betaindex[index], update.mu, beta,
                  lambda[index,,drop=FALSE], mu.block, sigma,
                  verbose=FALSE, compact=TRUE)
      lambda[index,] <- ss$lambda
      bound[index] <- ss$bound
      cache[[b]] <- list(beta=ss$beta[[1]], words=ss$beta.words, sigma=ss$sigma)

      updates <- updates + 1L
      rho <- (updates + control$tau)^(-control$kappa)
      scale <- N/length(index)
      S.beta <- (1 - rho) * S.beta
      S.beta[,ss$beta.words] <- S.beta[,ss$beta.words,drop=FALSE] +
        (rho * scale) * ss$beta[[1]]
      S.nu <- (1 - rho) * S.nu + (rho * scale) * ss$sigma
      S.lambda.sum <- (1 - rho) * S.lambda.sum +
        (rho * scale) * colSums(ss$lambda)
      S.lambda.cross <- (1 - rho) * S.lambda.cross +
        (rho * scale) * crossprod(ss$lambda)
      if(has.covariates) {
        S.x.lambda <- (1 - rho) * S.x.lambda +
          (rho * scale) * crossprod(X[index,,drop=FALSE], ss$lambda)
      }

      mapped <- incremental.mstep(
        S.beta + smoothing, S.nu, S.lambda.sum, S.lambda.cross,
        S.x.lambda, x.cross, N, settings, gamma.current=mu$gamma,
        update.gamma=TRUE
      )
      mu <- mapped$mu
      sigma <- mapped$sigma
      beta <- mapped$beta

      # Progress is measured by comparing each block's bound with its own
      # value on the previous pass: minibatch bounds vary far more across
      # blocks (which documents they hold) than across passes.
      per.token <- sum(ss$bound)/block.tokens[b]
      window <- c(window, per.token)
      if(!is.na(block.bound[b])) window.gain <- c(window.gain, per.token - block.bound[b])
      block.bound[b] <- per.token
      record$minibatch.bound <- c(record$minibatch.bound, per.token)
      record$minibatch.epoch <- c(record$minibatch.epoch,
                                  epoch - 1 + position/nblocks)

      if(length(window) >= checkpoint.every) {
        value <- mean(window)
        gain <- if(length(window.gain)) mean(window.gain) else NA_real_
        change <- gain/abs(value)
        window <- numeric()
        window.gain <- numeric()
        record$checkpoint.bound <- c(record$checkpoint.bound, value)
        record$checkpoint.epoch <- c(record$checkpoint.epoch,
                                     epoch - 1 + position/nblocks)
        record$checkpoint.change <- c(record$checkpoint.change, change)
        if(verbose) {
          cat(sprintf("  Stochastic pass %.2f: bound per token %.5f (step size %.3f%s)\n",
                      epoch - 1 + position/nblocks, value, rho,
                      if(is.na(change)) "" else sprintf(", gain since last pass %.2e", change)))
        }
        if(!is.na(change)) {
          ncheckpoints <- ncheckpoints + 1L
          flat <- if(change < control$tol) flat + 1L else 0L
          if(ncheckpoints >= control$min.checkpoints &&
             flat >= control$patience) {
            stop.now <- TRUE
            stop.reason <- "bound gain per pass fell below stochastic.tol"
          }
        }
      }
      if(stop.now) break
    }
    record$epochs <- epoch - 1 + position/nblocks
    if(stop.now) break
  }
  record$updates <- updates
  record$stop.reason <- stop.reason
  record$final.step.size <- (updates + control$tau)^(-control$kappa)
  if(verbose) {
    cat(sprintf("Stochastic phase finished after %.2f passes (%d updates): %s.\n",
                record$epochs, updates, stop.reason))
  }
  # Hand the last pass's block statistics to the memoized phase when they fit
  # in its memory budget; otherwise it refills a cache with its own blocks.
  cache.mb <- 8 * K * sum(vapply(cache, function(x) length(x$words), numeric(1)))/1024^2
  record$cache.mb <- cache.mb
  carry <- cache.mb <= settings$incremental$max.cache.mb
  stochastic.finish(documents, vocab, settings, mu, sigma, beta, lambda,
                    record, globaltime,
                    cache=if(carry) cache else NULL, blocks=blocks,
                    bound=bound)
}


# Either refine with the memoized incremental method or return the stochastic
# solution directly (with the usual optional theta pass).
stochastic.finish <- function(documents, vocab, settings, mu, sigma, beta,
                              lambda, record, globaltime, cache=NULL,
                              blocks=NULL, bound=NULL) {
  if(isTRUE(settings$stochastic$finish)) {
    return(incremental.control(
      documents, vocab, settings,
      start=list(mu=mu, sigma=sigma, beta=beta, lambda=lambda,
                 record=record, label="Stochastic", cache=cache,
                 blocks=blocks, bound=bound,
                 batch.size=settings$stochastic$batch.size),
      globaltime=globaltime
    ))
  }
  verbose <- settings$verbose
  X <- settings$covariates$X
  betaindex <- settings$covariates$betaindex
  convergence <- list(
    bound=numeric(), its=0L, epochs=record$epochs, updates=record$updates,
    converged=identical(record$stop.reason, "bound gain per pass fell below stochastic.tol"),
    stop.reason=record$stop.reason, theta.computed=settings$incremental$final.theta,
    stochastic=record
  )
  lambda.output <- NULL
  if(settings$incremental$final.theta) {
    if(verbose) cat("Performing final theta-only document pass...\n")
    mu.full <- if(is.null(mu$gamma)) mu$mu else t(as.matrix(X %*% mu$gamma))
    theta.time <- proc.time()
    final <- estep.lambda(documents, betaindex, !is.null(mu$gamma), beta,
                          lambda, mu.full, sigma, verbose, return.bound=TRUE)
    lambda.output <- final$lambda
    convergence$final.bound <- sum(final$bound)
    convergence$final.bound.per.token <-
      convergence$final.bound/sum(settings$dim$wcounts$x)
    convergence$theta.update.time <- unname((proc.time() - theta.time)[3])
    convergence$theta.updated <- TRUE
  } else {
    convergence$theta.updated <- FALSE
  }
  if(!is.null(mu$gamma)) mu$mu <- t(as.matrix(X %*% mu$gamma))
  out <- list(
    mu=mu, sigma=sigma, beta=list(logbeta=lapply(beta, safelog)),
    settings=settings, vocab=vocab, convergence=convergence,
    theta=if(is.null(lambda.output)) NULL else theta.from.lambda(lambda.output),
    eta=lambda.output, invsigma=solve(sigma),
    time=unname((proc.time() - globaltime)[3]),
    version=utils::packageDescription("stm")$Version,
    inference="Stochastic"
  )
  class(out) <- "STM"
  out
}
