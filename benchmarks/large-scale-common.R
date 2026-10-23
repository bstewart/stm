parse_benchmark_args <- function(defaults) {
  args <- commandArgs(trailingOnly=TRUE)
  out <- defaults
  if(!length(args)) return(out)
  if(any(!startsWith(args, "--")) || any(!grepl("=", args, fixed=TRUE))) {
    stop("Arguments must use --name=value syntax")
  }
  for(argument in args) {
    pieces <- strsplit(sub("^--", "", argument), "=", fixed=TRUE)[[1]]
    name <- pieces[1]
    if(!name %in% names(out)) stop("Unknown argument: --", name)
    value <- paste(pieces[-1], collapse="=")
    template <- out[[name]]
    out[[name]] <- if(is.logical(template)) {
      match.arg(tolower(value), c("true", "false")) == "true"
    } else if(is.integer(template)) {
      as.integer(value)
    } else if(is.numeric(template)) {
      as.numeric(value)
    } else {
      value
    }
  }
  out
}


softmax_rows <- function(eta) {
  eta <- cbind(eta, 0)
  eta <- eta - apply(eta, 1, max)
  theta <- exp(eta)
  theta/rowSums(theta)
}


hungarian_assignment <- function(cost) {
  cost <- as.matrix(cost)
  n <- nrow(cost)
  m <- ncol(cost)
  if(n > m) stop("Hungarian implementation requires rows <= columns")
  u <- numeric(n + 1L)
  v <- numeric(m + 1L)
  p <- integer(m + 1L)
  way <- integer(m + 1L)

  for(i in seq_len(n)) {
    p[1] <- i
    j0 <- 1L
    minv <- rep(Inf, m + 1L)
    used <- rep(FALSE, m + 1L)
    repeat {
      used[j0] <- TRUE
      i0 <- p[j0]
      delta <- Inf
      j1 <- 1L
      for(column in seq_len(m)) {
        j <- column + 1L
        if(!used[j]) {
          current <- cost[i0, column] - u[i0 + 1L] - v[j]
          if(current < minv[j]) {
            minv[j] <- current
            way[j] <- j0
          }
          if(minv[j] < delta) {
            delta <- minv[j]
            j1 <- j
          }
        }
      }
      for(j in seq_len(m + 1L)) {
        if(used[j]) {
          u[p[j] + 1L] <- u[p[j] + 1L] + delta
          v[j] <- v[j] - delta
        } else {
          minv[j] <- minv[j] - delta
        }
      }
      j0 <- j1
      if(p[j0] == 0L) break
    }
    repeat {
      j1 <- way[j0]
      p[j0] <- p[j1]
      j0 <- j1
      if(j0 == 1L) break
    }
  }

  assignment <- integer(n)
  for(column in seq_len(m)) {
    row <- p[column + 1L]
    if(row > 0L) assignment[row] <- column
  }
  assignment
}


topic_truth_metrics <- function(beta, theta, truth) {
  beta <- beta/rowSums(beta)
  affinity <- sqrt(beta) %*% t(sqrt(truth$beta))
  assignment <- hungarian_assignment(1 - affinity)
  matched.affinity <- affinity[cbind(seq_len(nrow(beta)), assignment)]
  hellinger <- sqrt(pmax(0, 1 - matched.affinity))

  anchor.recall <- vapply(seq_len(nrow(beta)), function(topic) {
    estimated.top <- order(beta[topic,], decreasing=TRUE)[
      seq_along(truth$anchors[[assignment[topic]]])
    ]
    mean(estimated.top %in% truth$anchors[[assignment[topic]]])
  }, numeric(1))

  theta.correlation <- vapply(seq_len(ncol(theta)), function(topic) {
    stats::cor(theta[,topic], truth$theta[,assignment[topic]])
  }, numeric(1))

  list(
    assignment=assignment,
    mean.hellinger=mean(hellinger),
    max.hellinger=max(hellinger),
    mean.anchor.recall=mean(anchor.recall),
    minimum.anchor.recall=min(anchor.recall),
    mean.theta.correlation=mean(theta.correlation),
    minimum.theta.correlation=min(theta.correlation)
  )
}


fresh_model_evaluation <- function(model, documents, truth) {
  beta <- lapply(model$beta$logbeta, exp)
  N <- length(documents)
  K <- model$settings$dim$K
  betaindex <- model$settings$covariates$betaindex
  X <- model$settings$covariates$X
  update.mu <- !is.null(model$mu$gamma)
  mu <- if(update.mu) {
    t(as.matrix(X %*% model$mu$gamma))
  } else {
    model$mu$mu
  }
  lambda <- model$eta
  if(is.null(lambda)) {
    lambda <- if(update.mu) {
      t(mu)
    } else {
      matrix(rep(as.numeric(mu), each=N), nrow=N, ncol=K - 1L)
    }
  }

  started <- proc.time()
  estep.arguments <- list(
    documents=documents, beta.index=betaindex, update.mu=update.mu,
    beta=beta, lambda.old=lambda, mu=mu, sigma=model$sigma,
    verbose=FALSE
  )
  if("compact" %in% names(formals(stm:::estep))) {
    estep.arguments$compact <- TRUE
  }
  evaluated <- do.call(stm:::estep, estep.arguments)
  elapsed <- unname((proc.time() - started)[3])
  theta <- softmax_rows(evaluated$lambda)
  total.tokens <- sum(vapply(documents, function(document) sum(document[2,]),
                             numeric(1)))
  truth.metrics <- topic_truth_metrics(beta[[1]], theta, truth)
  c(
    list(
      fresh.bound=sum(evaluated$bound),
      fresh.bound.per.token=sum(evaluated$bound)/total.tokens,
      evaluation.seconds=elapsed
    ),
    truth.metrics
  )
}
