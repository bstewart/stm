# Incremental (memoized) variational EM for finite corpora, in the sense of
# Neal and Hinton (1998).  Documents are split into fixed, shuffled blocks.
# A first pass fills a cache of per-block sufficient statistics; later
# block visits subtract the stale contribution, add the refreshed one and
# run the same M-step as batch inference.  There is no learning rate.

incremental.defaults <- function(N, K, V, max.cache.mb=256,
                                batch.size=NULL, max.blocks=12) {
  bytes.per.dense.block <- 8 * K * V
  memory.blocks <- floor(max.cache.mb * 1024^2 / bytes.per.dense.block)
  memory.blocks <- max(1, memory.blocks)

  if(is.null(batch.size)) {
    target.blocks <- min(max.blocks, N)
    nblocks <- min(target.blocks, memory.blocks, N)
    batch.size <- ceiling(N / nblocks)
    batch.size <- max(batch.size, min(N, max(128, 8 * K)))
  }
  batch.size <- min(N, as.integer(ceiling(batch.size)))
  nblocks <- ceiling(N / batch.size)

  max.epochs <- if(N < 1000) {
    60L
  } else if(N < 5000) {
    40L
  } else if(N < 20000) {
    30L
  } else {
    20L
  }

  eval.size <- if(N >= 10000) min(1000L, N) else 0L

  list(batch.size=batch.size, nblocks=nblocks,
       max.epochs=max.epochs, min.epochs=3L, patience=2L,
       final.theta=(N < 10000), max.cache.mb=max.cache.mb,
       progressive.fill=TRUE, eval.size=eval.size,
       eval.every=if(eval.size > 0) 1L else NULL,
       eval.patience=2L,
       early.stop=if(eval.size > 0) "both" else "cache",
       full.eval.every=NULL,
       estimated.cache.mb=nblocks * bytes.per.dense.block / 1024^2)
}


vb.variational.reg.stats <- function(XYcorr, YYcorr, Xcorr, N,
                                     b0=1, d0=1, maxits=1000) {
  Xcorr <- as.matrix(Xcorr)
  XYcorr <- as.numeric(XYcorr)
  an <- (1 + N)/2
  D <- ncol(Xcorr)
  w <- rep(0, D)
  error.prec <- 1
  converge <- Inf
  cn <- D
  dn <- 1
  Ea <- cn/dn
  ba <- 1
  ct <- 1L

  while(converge > .0001) {
    w.old <- w
    ppmat <- diag(c(0, rep(as.numeric(Ea), D - 1)), nrow=D)
    invV <- error.prec * Xcorr + ppmat
    V <- chol2inv(chol(invV))
    w <- as.numeric(error.prec * V %*% XYcorr)

    sse <- YYcorr - 2 * crossprod(w, XYcorr) +
      crossprod(w, Xcorr %*% w)
    sse <- max(0, as.numeric(sse))
    bn <- .5 * (sse + sum(diag(Xcorr %*% V))) + ba
    error.prec <- an/bn
    ba <- 1/(error.prec + b0)

    w0 <- w[1]
    w.penalized <- w[-1]
    da <- 2/(Ea + d0)
    dn <- 2 * da + crossprod(w.penalized) + sum(diag(V)[-1])
    Ea <- cn/dn
    w <- c(w0, w.penalized)
    ct <- ct + 1L
    if(ct > maxits) {
      stop("Prevalence regression failing to converge within iteration limit. ",
           "You can increase gamma.maxits in control.")
    }
    converge <- sum(abs(w - w.old))
  }
  w
}


incremental.mstep <- function(beta.counts, nu.sum, lambda.sum, lambda.cross,
                             x.lambda, x.cross, N, settings,
                             gamma.current=NULL, update.gamma=TRUE) {
  Kminus1 <- length(lambda.sum)
  mode <- settings$gamma$mode

  if(mode == "CTM") {
    mu.vector <- lambda.sum/N
    mu <- list(mu=matrix(mu.vector, ncol=1), gamma=NULL)
    covariance <- lambda.cross - N * tcrossprod(mu.vector)
  } else {
    gamma <- gamma.current
    if(update.gamma || is.null(gamma)) {
      gamma <- matrix(0, nrow=nrow(x.lambda), ncol=Kminus1)
      for(k in seq_len(Kminus1)) {
        gamma[,k] <- vb.variational.reg.stats(
          XYcorr=x.lambda[,k], YYcorr=lambda.cross[k,k],
          Xcorr=x.cross, N=N, maxits=settings$gamma$maxits
        )
      }
    }
    covariance <- lambda.cross - crossprod(gamma, x.lambda) -
      crossprod(x.lambda, gamma) + crossprod(gamma, x.cross %*% gamma)
    mu <- list(mu=NULL, gamma=gamma)
  }

  sigma <- (covariance + nu.sum)/N
  sigma <- (sigma + t(sigma))/2
  sigma <- diag(diag(sigma), nrow=Kminus1) * settings$sigma$prior +
    (1 - settings$sigma$prior) * sigma
  if(inherits(try(chol.default(sigma), silent=TRUE), "try-error")) {
    eig <- eigen(sigma, symmetric=TRUE)
    floor.value <- max(1e-8, max(eig$values) * 1e-10)
    eig$values[eig$values < floor.value] <- floor.value
    sigma <- tcrossprod(sweep(eig$vectors, 2, sqrt(eig$values), `*`))
    sigma <- (sigma + t(sigma))/2
  }

  beta.counts[beta.counts < 0 & beta.counts > -1e-10] <- 0
  if(any(beta.counts < 0)) {
    stop("Negative memoized topic-word sufficient statistic encountered")
  }
  totals <- rowSums(beta.counts)
  if(any(!is.finite(totals)) || any(totals <= 0)) {
    stop("Invalid memoized topic-word sufficient statistics")
  }
  beta <- list(beta.counts/totals)

  list(mu=mu, sigma=sigma, beta=beta)
}


incremental.sample <- function(n, size=n, seed) {
  has.seed <- exists(".Random.seed", envir=.GlobalEnv, inherits=FALSE)
  if(has.seed) old.seed <- get(".Random.seed", envir=.GlobalEnv,
                               inherits=FALSE)
  on.exit({
    if(has.seed) {
      assign(".Random.seed", old.seed, envir=.GlobalEnv)
    } else if(exists(".Random.seed", envir=.GlobalEnv, inherits=FALSE)) {
      rm(".Random.seed", envir=.GlobalEnv)
    }
  })
  set.seed(seed)
  sample.int(n, size=size, replace=FALSE)
}


incremental.seed.offset <- function(seed, offset) {
  as.integer((as.double(seed) + as.double(offset)) %%
               .Machine$integer.max)
}


incremental.blocks <- function(N, batch.size, seed) {
  ordering <- incremental.sample(N, seed=seed)
  split(ordering, ceiling(seq_along(ordering)/batch.size))
}


incremental.epoch.order <- function(nblocks, seed, epoch) {
  epoch.seed <- incremental.seed.offset(seed, 1009 * as.double(epoch))
  incremental.sample(nblocks, seed=epoch.seed)
}


theta.from.lambda <- function(lambda) {
  eta <- cbind(lambda, 0)
  eta <- eta - matrixStats::rowMaxs(eta)
  theta <- exp(eta)
  theta/rowSums(theta)
}


incremental.evaluate.bound <- function(documents, index, betaindex, X, mu,
                                      beta, lambda, sigma) {
  update.mu <- !is.null(mu$gamma)
  mu.eval <- if(update.mu) {
    t(as.matrix(X[index,,drop=FALSE] %*% mu$gamma))
  } else {
    mu$mu
  }
  evaluated <- estep(
    documents[index], betaindex[index], update.mu, beta,
    lambda[index,,drop=FALSE], mu.eval, sigma, verbose=FALSE, compact=TRUE
  )
  tokens <- sum(vapply(documents[index], function(document) sum(document[2,]),
                       numeric(1)))
  total <- sum(evaluated$bound)
  list(bound=total, per.token=total/tokens)
}


incremental.relative.change <- function(current, previous) {
  if(length(previous) == 0 || !is.finite(current) || !is.finite(previous)) {
    return(NA_real_)
  }
  (current - previous)/abs(previous)
}


# start: optional list(mu, sigma, beta, lambda, record, label) supplying the
# starting values (used by stochastic.control to hand over to the memoized
# phase); globaltime: optional earlier start time so the reported time covers
# both phases.
incremental.control <- function(documents, vocab, settings, model=NULL,
                                start=NULL, globaltime=proc.time()) {
  N <- settings$dim$N
  K <- settings$dim$K
  V <- settings$dim$V
  verbose <- settings$verbose
  betaindex <- settings$covariates$betaindex
  X <- settings$covariates$X
  control <- settings$incremental
  label <- if(is.null(start$label)) "Incremental" else start$label

  prior.epochs <- 0L
  prior.updates <- 0L
  prior.convergence <- NULL
  if(!is.null(start)) {
    mu <- start$mu
    sigma <- start$sigma
    beta <- start$beta
    lambda <- start$lambda
  } else if(is.null(model)) {
    if(verbose) cat(switch(EXPR=settings$init$mode,
      Spectral="Beginning Spectral Initialization \n",
      LDA="Beginning LDA Initialization \n",
      Random="Beginning Random Initialization \n",
      Custom="Beginning Custom Initialization \n"))
    initialized <- stm.init(documents, settings)
    if(settings$dim$K == 0) {
      stop("K=0 topic selection is not supported with incremental inference")
    }
    mu <- list(mu=initialized$mu, gamma=NULL)
    sigma <- initialized$sigma
    beta <- initialized$beta
    lambda <- initialized$lambda
  } else {
    if(!inherits(model, "STM")) stop("model must inherit from class 'STM'")
    if(!identical(model$inference, "Incremental")) {
      stop("Incremental inference can only restart a incremental STM fit")
    }
    old.dim <- model$settings$dim
    if(old.dim$N != N || old.dim$K != K || old.dim$V != V) {
      stop("Restarted incremental model dimensions do not match the supplied corpus")
    }
    if(verbose) cat("Restarting with incremental inference...\n")
    mu <- model$mu
    sigma <- model$sigma
    beta <- lapply(model$beta$logbeta, exp)
    lambda <- model$eta
    if(is.null(lambda)) {
      if(is.null(mu$gamma)) {
        lambda <- matrix(rep(as.numeric(mu$mu), each=N), nrow=N)
      } else {
        lambda <- as.matrix(X %*% mu$gamma)
      }
    }
    prior.convergence <- model$convergence
    if(!is.null(prior.convergence$its)) {
      prior.epochs <- as.integer(prior.convergence$its)
    }
    if(!is.null(prior.convergence$updates)) {
      prior.updates <- as.integer(prior.convergence$updates)
    }
    if(control$max.epochs <= prior.epochs) {
      stop("incremental.max.epochs is a total epoch limit when restarting and must exceed the previous run")
    }
  }

  epochs.to.run <- control$max.epochs - prior.epochs
  if(control$max.epochs == 0) {
    convergence <- list(
      bound=numeric(), bound.epoch=numeric(), relative.change=numeric(),
      its=0L, epochs=0, updates=0L, converged=FALSE,
      stop.reason="zero epochs requested", theta.computed=TRUE
    )
    lambda.output <- lambda
  } else {
    carried <- !is.null(start$cache)
    if(carried) {
      blocks <- start$blocks
      control$batch.size <- start$batch.size
      settings$incremental$batch.size <- start$batch.size
    } else {
      blocks <- incremental.blocks(N, control$batch.size, control$seed)
    }
    nblocks <- length(blocks)
    eval.index <- if(control$eval.size > 0) {
      incremental.sample(N, size=control$eval.size, seed=control$eval.seed)
    } else {
      integer()
    }
    settings$incremental$eval.indices <- eval.index

    cache <- vector("list", nblocks)
    beta.counts <- matrix(0, nrow=K, ncol=V)
    nu.sum <- matrix(0, nrow=K-1, ncol=K-1)
    lambda.sum <- numeric(K-1)
    lambda.cross <- matrix(0, nrow=K-1, ncol=K-1)
    x.lambda <- if(is.null(X)) NULL else matrix(0, nrow=ncol(X), ncol=K-1)
    x.cross <- if(is.null(X)) NULL else crossprod(X)
    bound <- rep(NA_real_, N)

    history <- if(is.null(prior.convergence$bound)) numeric() else
      prior.convergence$bound
    bound.epoch <- if(is.null(prior.convergence$bound.epoch)) {
      seq_along(history)
    } else prior.convergence$bound.epoch
    relative.change <- if(is.null(prior.convergence$relative.change)) {
      numeric()
    } else prior.convergence$relative.change
    eval.bound <- if(is.null(prior.convergence$eval.bound)) numeric() else
      prior.convergence$eval.bound
    eval.per.token <- if(is.null(prior.convergence$eval.per.token)) numeric() else
      prior.convergence$eval.per.token
    eval.epoch <- if(is.null(prior.convergence$eval.epoch)) numeric() else
      prior.convergence$eval.epoch
    eval.relative.change <- if(is.null(prior.convergence$eval.relative.change)) {
      numeric()
    } else prior.convergence$eval.relative.change
    full.bound <- if(is.null(prior.convergence$full.bound)) numeric() else
      prior.convergence$full.bound
    full.per.token <- if(is.null(prior.convergence$full.per.token)) numeric() else
      prior.convergence$full.per.token
    full.epoch <- if(is.null(prior.convergence$full.epoch)) numeric() else
      prior.convergence$full.epoch

    cache.patience.count <- 0L
    eval.patience.count <- 0L
    converged <- FALSE
    stop.reason <- "maximum epochs reached"
    updates <- prior.updates
    run.updates <- 0L
    epoch.equivalent <- as.numeric(prior.epochs)

    record.monitoring <- function(epoch.number, epoch.value, force=FALSE) {
      evaluate.now <- length(eval.index) > 0 &&
        (force || epoch.number %% control$eval.every == 0L)
      if(evaluate.now) {
        evaluated <- incremental.evaluate.bound(
          documents, eval.index, betaindex, X, mu, beta, lambda, sigma
        )
        previous <- if(length(eval.per.token)) tail(eval.per.token, 1) else
          numeric()
        eval.bound <<- c(eval.bound, evaluated$bound)
        eval.per.token <<- c(eval.per.token, evaluated$per.token)
        eval.epoch <<- c(eval.epoch, epoch.value)
        change <- incremental.relative.change(evaluated$per.token, previous)
        eval.relative.change <<- c(eval.relative.change, change)
        stable <- is.finite(change) &&
          abs(change) < settings$convergence$em.converge.thresh
        eval.patience.count <<- if(stable) eval.patience.count + 1L else 0L
        if(verbose) {
          cat(sprintf(paste0("  Fresh evaluation bound = %.6f per token ",
                             "(relative change = %.3g).\n"),
                      evaluated$per.token, change))
        }
      }

      full.now <- !is.null(control$full.eval.every) &&
        (force || epoch.number %% control$full.eval.every == 0L)
      if(full.now) {
        evaluated <- incremental.evaluate.bound(
          documents, seq_len(N), betaindex, X, mu, beta, lambda, sigma
        )
        full.bound <<- c(full.bound, evaluated$bound)
        full.per.token <<- c(full.per.token, evaluated$per.token)
        full.epoch <<- c(full.epoch, epoch.value)
        if(verbose) {
          cat(sprintf("  Fresh full-corpus bound = %.6f per token.\n",
                      evaluated$per.token))
        }
      }
      invisible(NULL)
    }

    check.convergence <- function() {
      old.enough <- epoch.equivalent >= control$min.epochs
      cache.stop <- old.enough &&
        cache.patience.count >= control$patience
      eval.stop <- old.enough && length(eval.index) > 0 &&
        eval.patience.count >= control$eval.patience
      switch(control$early.stop,
             cache=cache.stop,
             evaluation=eval.stop,
             both=cache.stop && eval.stop)
    }

    if(verbose) {
      cat(sprintf("Beginning %s inference with %d blocks of about %d documents.\n",
                  if(is.null(start)) "incremental" else "memoized refinement",
                  nblocks, control$batch.size))
      if(!carried) cat("Populating the memoized sufficient-statistic cache...\n")
    }

    if(carried) {
    # The stochastic phase visited every block during its last pass and kept
    # each block's statistics, so the cache is already complete: sum it and
    # take one exact M-step instead of spending a pass on refilling it.
    cache <- start$cache
    for(b in seq_along(blocks)) {
      beta.counts[,cache[[b]]$words] <-
        beta.counts[,cache[[b]]$words,drop=FALSE] + cache[[b]]$beta
      nu.sum <- nu.sum + cache[[b]]$sigma
    }
    lambda.sum <- colSums(lambda)
    lambda.cross <- crossprod(lambda)
    if(!is.null(X)) x.lambda <- crossprod(X, lambda)
    bound <- start$bound
    } else {
    # The first pass fills the cache with an exact E-step for every block.
    # Unless this is a restart of a fitted model, the global parameters are
    # also updated after each block: blocks not yet visited contribute the
    # statistics implied by the starting values (their token mass spread over
    # the starting beta, the starting covariance for each document and their
    # stored document means), so the M-step is well defined from the first
    # block on and the statistics become exact as the pass completes.
    progressive <- is.null(model) && !carried && isTRUE(control$progressive.fill)
    if(progressive) {
      ntokens <- sum(settings$dim$wcounts$x)
      block.tokens <- vapply(blocks, function(index) {
        sum(vapply(documents[index], function(document) sum(document[2,]),
                   numeric(1)))
      }, numeric(1))
      beta.start <- beta[[1]]
      sigma.start <- sigma
      unseen.tokens <- ntokens
      unseen.docs <- N
      lambda.sum <- colSums(lambda)
      lambda.cross <- crossprod(lambda)
      if(!is.null(X)) x.lambda <- crossprod(X, lambda)
    }
    for(b in seq_along(blocks)) {
      index <- blocks[[b]]
      update.mu <- !is.null(mu$gamma)
      mu.block <- if(update.mu) {
        t(as.matrix(X[index,,drop=FALSE] %*% mu$gamma))
      } else {
        mu$mu
      }
      old.lambda <- lambda[index,,drop=FALSE]
      ss <- estep(documents[index], betaindex[index], update.mu,
                  beta, old.lambda, mu.block, sigma,
                  verbose=FALSE, compact=TRUE)

      beta.counts[,ss$beta.words] <-
        beta.counts[,ss$beta.words,drop=FALSE] + ss$beta[[1]]
      nu.sum <- nu.sum + ss$sigma
      lambda[index,] <- ss$lambda
      if(progressive) {
        lambda.sum <- lambda.sum - colSums(old.lambda) + colSums(ss$lambda)
        lambda.cross <- lambda.cross - crossprod(old.lambda) +
          crossprod(ss$lambda)
        if(!is.null(X)) {
          X.block <- X[index,,drop=FALSE]
          x.lambda <- x.lambda - crossprod(X.block, old.lambda) +
            crossprod(X.block, ss$lambda)
        }
      } else {
        lambda.sum <- lambda.sum + colSums(ss$lambda)
        lambda.cross <- lambda.cross + crossprod(ss$lambda)
        if(!is.null(X)) {
          x.lambda <- x.lambda + crossprod(X[index,,drop=FALSE], ss$lambda)
        }
      }
      bound[index] <- ss$bound
      cache[[b]] <- list(beta=ss$beta[[1]], words=ss$beta.words,
                         sigma=ss$sigma)

      if(progressive && b < nblocks) {
        unseen.tokens <- unseen.tokens - block.tokens[b]
        unseen.docs <- unseen.docs - length(index)
        mapped <- incremental.mstep(
          beta.counts + beta.start * (unseen.tokens/K),
          nu.sum + unseen.docs * sigma.start,
          lambda.sum, lambda.cross, x.lambda, x.cross,
          N, settings, gamma.current=mu$gamma, update.gamma=TRUE
        )
        mu <- mapped$mu
        sigma <- mapped$sigma
        beta <- mapped$beta
      }
    }
    }

    mapped <- incremental.mstep(
      beta.counts, nu.sum, lambda.sum, lambda.cross, x.lambda, x.cross,
      N, settings, gamma.current=mu$gamma, update.gamma=TRUE
    )
    mu <- mapped$mu
    sigma <- mapped$sigma
    beta <- mapped$beta
    current.bound <- sum(bound)
    previous.bound <- if(length(history)) tail(history, 1) else numeric()
    history <- c(history, current.bound)
    epoch.equivalent <- prior.epochs + 1
    bound.epoch <- c(bound.epoch, epoch.equivalent)
    change <- incremental.relative.change(current.bound, previous.bound)
    relative.change <- c(relative.change, change)
    stable <- is.finite(change) &&
      abs(change) < settings$convergence$em.converge.thresh
    cache.patience.count <- if(stable) cache.patience.count + 1L else 0L
    record.monitoring(prior.epochs + 1L, epoch.equivalent)

    if(verbose) {
      cat(sprintf(paste0("Completed incremental epoch %d ",
                         "(approximate bound = %.3f, relative change = %.3g).\n"),
                  prior.epochs + 1L, current.bound, change))
    }
    converged <- check.convergence()

    if(epochs.to.run > 1 && !converged) {
      for(run.epoch in 2:epochs.to.run) {
        absolute.epoch <- prior.epochs + run.epoch
        block.order <- incremental.epoch.order(nblocks, control$seed,
                                               absolute.epoch)
        for(position in seq_along(block.order)) {
          b <- block.order[position]
          index <- blocks[[b]]
          old.lambda <- lambda[index,,drop=FALSE]
          update.mu <- !is.null(mu$gamma)
          mu.block <- if(update.mu) {
            t(as.matrix(X[index,,drop=FALSE] %*% mu$gamma))
          } else {
            mu$mu
          }
          ss <- estep(documents[index], betaindex[index], update.mu,
                      beta, old.lambda, mu.block, sigma,
                      verbose=FALSE, compact=TRUE)

          old <- cache[[b]]
          if(!identical(old$words, ss$beta.words)) {
            stop("Internal incremental block vocabulary changed")
          }
          beta.counts[,ss$beta.words] <-
            beta.counts[,ss$beta.words,drop=FALSE] - old$beta + ss$beta[[1]]
          nu.sum <- nu.sum - old$sigma + ss$sigma
          lambda.sum <- lambda.sum - colSums(old.lambda) + colSums(ss$lambda)
          lambda.cross <- lambda.cross - crossprod(old.lambda) +
            crossprod(ss$lambda)
          if(!is.null(X)) {
            X.block <- X[index,,drop=FALSE]
            x.lambda <- x.lambda - crossprod(X.block, old.lambda) +
              crossprod(X.block, ss$lambda)
          }
          lambda[index,] <- ss$lambda
          bound[index] <- ss$bound
          cache[[b]] <- list(beta=ss$beta[[1]], words=ss$beta.words,
                             sigma=ss$sigma)
          updates <- updates + 1L
          run.updates <- run.updates + 1L

          update.gamma <- is.null(X) ||
            updates %% control$gamma.update.every == 0L
          mapped <- incremental.mstep(
            beta.counts, nu.sum, lambda.sum, lambda.cross, x.lambda, x.cross,
            N, settings, gamma.current=mu$gamma,
            update.gamma=update.gamma
          )
          mu <- mapped$mu
          sigma <- mapped$sigma
          beta <- mapped$beta
        }

        current.bound <- sum(bound)
        previous.bound <- tail(history, 1)
        history <- c(history, current.bound)
        epoch.equivalent <- prior.epochs + 1 + run.updates/nblocks
        bound.epoch <- c(bound.epoch, epoch.equivalent)
        change <- incremental.relative.change(current.bound, previous.bound)
        relative.change <- c(relative.change, change)
        stable <- is.finite(change) &&
          abs(change) < settings$convergence$em.converge.thresh
        cache.patience.count <- if(stable) cache.patience.count + 1L else 0L
        record.monitoring(absolute.epoch, epoch.equivalent)

        if(verbose && (run.epoch == 2L ||
                       absolute.epoch %% settings$topicreportevery == 0L)) {
          cat(sprintf(paste0("Completed incremental epoch %.3f ",
                             "(approximate bound = %.3f, relative change = %.3g).\n"),
                      epoch.equivalent, current.bound, change))
        }
        converged <- check.convergence()
        if(converged) {
          stop.reason <- paste0(control$early.stop,
                                " bound convergence criterion reached")
          break
        }
      }
    } else if(converged) {
      stop.reason <- paste0(control$early.stop,
                            " bound convergence criterion reached")
    }

    # Ensure periodically updated prevalence parameters correspond to the
    # final cached moments even when a custom update interval does not divide
    # the number of completed block updates.
    mapped <- incremental.mstep(
      beta.counts, nu.sum, lambda.sum, lambda.cross, x.lambda, x.cross,
      N, settings, gamma.current=mu$gamma, update.gamma=TRUE
    )
    mu <- mapped$mu
    sigma <- mapped$sigma
    beta <- mapped$beta
    convergence <- list(
      bound=history, bound.epoch=bound.epoch,
      relative.change=relative.change,
      eval.bound=eval.bound, eval.per.token=eval.per.token,
      eval.epoch=eval.epoch, eval.relative.change=eval.relative.change,
      full.bound=full.bound, full.per.token=full.per.token,
      full.epoch=full.epoch,
      its=length(history), epochs=epoch.equivalent,
      updates=updates, converged=converged, stop.reason=stop.reason,
      theta.computed=control$final.theta,
      blocks=nblocks, batch.size=control$batch.size,
      cache.mb=as.numeric(utils::object.size(cache))/1024^2,
      cache.word.columns=sum(vapply(cache, function(x) length(x$words),
                                    numeric(1)))
    )

    if(control$final.theta) {
      if(verbose) cat("Performing final theta-only document pass...\n")
      mu.full <- if(is.null(mu$gamma)) mu$mu else
        t(as.matrix(X %*% mu$gamma))
      theta.time <- proc.time()
      final <- estep.lambda(
        documents, betaindex, !is.null(mu$gamma), beta, lambda,
        mu.full, sigma, verbose, return.bound=TRUE
      )
      lambda.output <- final$lambda
      convergence$final.bound <- sum(final$bound)
      convergence$final.bound.per.token <-
        convergence$final.bound/sum(settings$dim$wcounts$x)
      convergence$theta.update.time <-
        unname((proc.time() - theta.time)[3])
      convergence$theta.updated <- TRUE
    } else {
      lambda.output <- NULL
      convergence$theta.updated <- FALSE
    }
  }

  if(!is.null(mu$gamma)) {
    mu$mu <- t(as.matrix(X %*% mu$gamma))
  }
  beta.out <- list(logbeta=lapply(beta, safelog))
  if(!is.null(lambda.output)) {
    theta <- theta.from.lambda(lambda.output)
  } else {
    theta <- NULL
  }
  if(!is.null(start$record)) convergence$stochastic <- start$record
  elapsed <- unname((proc.time() - globaltime)[3])
  out <- list(
    mu=mu, sigma=sigma, beta=beta.out, settings=settings,
    vocab=vocab, convergence=convergence, theta=theta,
    eta=lambda.output, invsigma=solve(sigma), time=elapsed,
    version=utils::packageDescription("stm")$Version,
    inference=label
  )
  class(out) <- "STM"
  out
}


#' Compute document-topic proportions for a fitted model
#'
#' Runs the document-level variational update for every document at the
#' model's current (fixed) global parameters and stores the resulting
#' \code{theta} and \code{eta} in the model.  Topic-word distributions, the
#' prevalence regression and the topic covariance are not changed.
#'
#' This is mainly for models fit with \code{inference="Incremental"} on large
#' corpora (10,000 or more documents), where the final document pass is
#' skipped by default and \code{theta} is \code{NULL} until this function is
#' run.  Functions that need \code{theta}, such as \code{\link{estimateEffect}}
#' and \code{\link{findThoughts}}, will say so.  It can also be applied to
#' any \code{STM} object to recompute \code{theta} with the stored global
#' parameters; note that batch models store the \code{theta} from their last
#' E-step, which precedes the final M-step, so the recomputed values will
#' differ somewhat.  Unlike \code{\link{fitNewDocuments}}, which handles new
#' documents and new covariate values, this function uses the prevalence
#' covariates stored in the model and requires the original training documents.
#'
#' @param model An \code{STM} object.
#' @param documents The original corpus in the native \pkg{stm} list format
#' (the same \code{documents} object the model was fit to).
#' @param verbose Logical; print progress of the document pass.
#' @return The model with \code{theta} and \code{eta} filled in. The
#' \code{convergence} element gains \code{final.bound} (the approximate
#' variational bound at the fixed global parameters),
#' \code{final.bound.per.token}, \code{theta.updated} and
#' \code{theta.update.time}.
#' @seealso \code{\link{stm}}, \code{\link{fitNewDocuments}}
#' @examples
#' \donttest{
#' temp<-textProcessor(documents=gadarian$open.ended.response,metadata=gadarian)
#' out <- prepDocuments(temp$documents, temp$vocab, temp$meta)
#' mod <- stm(out$documents, out$vocab, K=3, inference="Incremental",
#'            max.em.its=5, control=list(incremental.final.theta=FALSE))
#' is.null(mod$theta)
#' mod <- updateTheta(mod, out$documents)
#' dim(mod$theta)
#' }
#' @export
updateTheta <- function(model, documents, verbose=TRUE) {
  if(!inherits(model, "STM")) stop("model must inherit from class 'STM'")
  if(!is.list(documents) || length(documents) != model$settings$dim$N) {
    stop("documents must be the original STM-format document list")
  }
  beta <- lapply(model$beta$logbeta, exp)
  update.mu <- !is.null(model$mu$gamma)
  mu <- model$mu$mu
  lambda.old <- model$eta
  if(is.null(lambda.old)) {
    if(update.mu) {
      lambda.old <- t(mu)
    } else {
      lambda.old <- matrix(rep(as.numeric(mu), each=length(documents)),
                           nrow=length(documents))
    }
  }
  started <- proc.time()
  refreshed <- estep.lambda(
    documents, model$settings$covariates$betaindex, update.mu,
    beta, lambda.old, mu, model$sigma, verbose, return.bound=TRUE
  )
  lambda <- refreshed$lambda
  model$eta <- lambda
  model$theta <- theta.from.lambda(lambda)
  model$convergence$theta.computed <- TRUE
  model$convergence$theta.updated <- TRUE
  model$convergence$final.bound <- sum(refreshed$bound)
  model$convergence$final.bound.per.token <-
    model$convergence$final.bound/sum(model$settings$dim$wcounts$x)
  model$convergence$theta.update.time <-
    unname((proc.time() - started)[3])
  model
}
