#!/usr/bin/env Rscript

# Component and policy ablations for stochastic STM inference.

suppressPackageStartupMessages(library(stm))
data(poliblog5k, package="stm")

out <- prepDocuments(
  poliblog5k.docs[seq_len(1000)], poliblog5k.voc,
  poliblog5k.meta[seq_len(1000),,drop=FALSE],
  lower.thresh=0, verbose=FALSE
)

# Compact versus dense block sufficient statistics at fixed globals.
initial <- stm(
  out$documents, out$vocab, K=20, init.type="Spectral", seed=91,
  max.em.its=0, verbose=FALSE
)
beta <- lapply(initial$beta$logbeta, exp)
index <- seq_len(160)
time_estep <- function(compact, repetitions=5L) {
  times <- numeric(repetitions)
  result <- NULL
  for(i in seq_len(repetitions)) {
    times[i] <- system.time(result <- stm:::estep(
      out$documents[index], rep(1, length(index)), FALSE, beta,
      initial$eta[index,,drop=FALSE], initial$mu$mu, initial$sigma,
      FALSE, compact=compact
    ))[["elapsed"]]
  }
  list(time=median(times), result=result)
}
dense.run <- time_estep(FALSE)
compact.run <- time_estep(TRUE)
dense.time <- dense.run$time
compact.time <- compact.run$time
dense <- dense.run$result
compact <- compact.run$result
compact.result <- data.frame(
  component="block_statistics_poliblog", dense_seconds=dense.time,
  optimized_seconds=compact.time, speedup=dense.time/compact.time,
  dense_beta_mb=as.numeric(object.size(dense$beta))/1024^2,
  compact_beta_mb=as.numeric(object.size(compact$beta))/1024^2,
  max_difference=max(abs(dense$beta[[1]][,compact$beta.words] -
                           compact$beta[[1]]))
)

# The poliblog vocabulary is small and a 160-document block covers most words.
# Repeat the allocation comparison in a large-vocabulary sparse block.
set.seed(99)
large.D <- 160L
large.V <- 50000L
large.K <- 20L
large.beta <- matrix(rexp(large.K * large.V), nrow=large.K)
large.beta <- large.beta/rowSums(large.beta)
large.documents <- lapply(seq_len(large.D), function(i) {
  rbind(sort(sample.int(large.V, 100)), sample.int(5, 100, TRUE))
})
large.lambda <- matrix(0, nrow=large.D, ncol=large.K-1)
time_large <- function(compact, repetitions=5L) {
  times <- numeric(repetitions)
  result <- NULL
  for(i in seq_len(repetitions)) {
    times[i] <- system.time(result <- stm:::estep(
      large.documents, rep(1, large.D), FALSE, list(large.beta),
      large.lambda, matrix(0, large.K-1, 1), diag(2, large.K-1),
      FALSE, compact=compact
    ))[["elapsed"]]
  }
  list(time=median(times), result=result)
}
large.dense <- time_large(FALSE)
large.compact <- time_large(TRUE)
large.result <- data.frame(
  component="block_statistics_large_vocab",
  dense_seconds=large.dense$time,
  optimized_seconds=large.compact$time,
  speedup=large.dense$time/large.compact$time,
  dense_beta_mb=as.numeric(object.size(large.dense$result$beta))/1024^2,
  compact_beta_mb=as.numeric(object.size(large.compact$result$beta))/1024^2,
  max_difference=max(abs(
    large.dense$result$beta[[1]][,large.compact$result$beta.words] -
      large.compact$result$beta[[1]]
  ))
)

fresh_bound <- function(model, update.mu=FALSE) {
  sum(stm:::estep(
    out$documents, model$settings$covariates$betaindex, update.mu,
    lapply(model$beta$logbeta, exp), model$eta,
    model$mu$mu, model$sigma, FALSE
  )$bound)
}

# End-to-end block-size policy at a fixed seed.
block.results <- lapply(c(40L, 80L, 160L, 320L), function(batch.size) {
  elapsed <- system.time(model <- stm(
    out$documents, out$vocab, K=20, init.type="Spectral", seed=91,
    inference="Incremental", emtol=1e-5, verbose=FALSE,
    control=list(incremental.batch.size=batch.size)
  ))[["elapsed"]]
  data.frame(
    batch_size=batch.size, elapsed=unname(elapsed),
    epochs=model$convergence$its, updates=model$convergence$updates,
    converged=model$convergence$converged, bound=fresh_bound(model),
    cache_mb=model$convergence$cache.mb
  )
})
block.results <- do.call(rbind, block.results)

# Existing dense ngroups comparator at approximately the default block count.
legacy.time <- system.time(legacy <- stm(
  out$documents, out$vocab, K=20, init.type="Spectral", seed=91,
  ngroups=7, max.em.its=120, emtol=1e-5, verbose=FALSE
))[["elapsed"]]
legacy.result <- data.frame(
  method=c("legacy_ngroups", "new_stochastic"),
  elapsed=c(unname(legacy.time),
            block.results$elapsed[block.results$batch_size == 160]),
  epochs=c(legacy$convergence$its,
           block.results$epochs[block.results$batch_size == 160]),
  bound=c(fresh_bound(legacy),
          block.results$bound[block.results$batch_size == 160])
)

# Pooled prevalence update frequency: every block versus once per epoch.
prevalence.results <- lapply(c(1L, 8L), function(update.every) {
  elapsed <- system.time(model <- stm(
    out$documents, out$vocab, K=15,
    prevalence=~rating + blog, data=out$meta,
    init.type="Spectral", seed=92, inference="Incremental",
    emtol=1e-5, verbose=FALSE,
    control=list(incremental.batch.size=128,
                 incremental.gamma.update.every=update.every)
  ))[["elapsed"]]
  data.frame(
    gamma_update_every=update.every, elapsed=unname(elapsed),
    epochs=model$convergence$its, bound=fresh_bound(model, TRUE)
  )
})
prevalence.results <- do.call(rbind, prevalence.results)

# Theta-only final pass versus a full E-step at identical globals.
model <- stm(
  out$documents, out$vocab, K=20, init.type="Spectral", seed=91,
  inference="Incremental", emtol=1e-5, verbose=FALSE
)
beta <- lapply(model$beta$logbeta, exp)
time_final <- function(theta.only, repetitions=5L) {
  times <- numeric(repetitions)
  result <- NULL
  for(i in seq_len(repetitions)) {
    times[i] <- system.time(result <- if(theta.only) {
      stm:::estep.lambda(
        out$documents, model$settings$covariates$betaindex, FALSE, beta,
        model$eta, model$mu$mu, model$sigma, FALSE
      )
    } else {
      stm:::estep(
        out$documents, model$settings$covariates$betaindex, FALSE, beta,
        model$eta, model$mu$mu, model$sigma, FALSE
      )
    })[["elapsed"]]
  }
  list(time=median(times), result=result)
}
lambda.run <- time_final(TRUE)
full.run <- time_final(FALSE)
lambda.time <- lambda.run$time
full.time <- full.run$time
lambda <- lambda.run$result
full <- full.run$result
theta.result <- data.frame(
  component="final_theta", dense_seconds=full.time,
  optimized_seconds=lambda.time, speedup=full.time/lambda.time,
  dense_beta_mb=NA_real_, compact_beta_mb=NA_real_,
  max_difference=max(abs(lambda - full$lambda))
)

print(rbind(compact.result, large.result, theta.result), row.names=FALSE)
print(block.results, row.names=FALSE)
print(prevalence.results, row.names=FALSE)
print(legacy.result, row.names=FALSE)
