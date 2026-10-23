make_incremental_case <- function(n=80, K=6, seed=17) {
  prepared <- prepDocuments(
    poliblog5k.docs[seq_len(n)], poliblog5k.voc,
    poliblog5k.meta[seq_len(n),,drop=FALSE],
    lower.thresh=0, verbose=FALSE
  )
  set.seed(seed)
  beta <- matrix(rgamma(K * length(prepared$vocab), .2), nrow=K)
  beta <- beta/rowSums(beta)
  list(
    documents=prepared$documents, vocab=prepared$vocab, meta=prepared$meta,
    beta=list(beta), lambda=matrix(0, nrow=n, ncol=K-1),
    mu=matrix(0, nrow=K-1, ncol=1), sigma=diag(2, K-1)
  )
}

test_that("compact E-step statistics reproduce dense statistics", {
  case <- make_incremental_case()
  index <- seq_len(12)
  dense <- stm:::estep(
    case$documents[index], rep(1, length(index)), FALSE,
    case$beta, case$lambda[index,,drop=FALSE], case$mu, case$sigma,
    FALSE, compact=FALSE
  )
  compact <- stm:::estep(
    case$documents[index], rep(1, length(index)), FALSE,
    case$beta, case$lambda[index,,drop=FALSE], case$mu, case$sigma,
    FALSE, compact=TRUE
  )
  restored <- matrix(0, nrow=nrow(case$beta[[1]]),
                     ncol=ncol(case$beta[[1]]))
  restored[,compact$beta.words] <- compact$beta[[1]]

  expect_equal(compact$lambda, dense$lambda, tolerance=1e-12)
  expect_equal(compact$sigma, dense$sigma, tolerance=1e-12)
  expect_equal(compact$bound, dense$bound, tolerance=1e-12)
  expect_equal(restored, dense$beta[[1]], tolerance=1e-12)
  expect_lt(ncol(compact$beta[[1]]), ncol(dense$beta[[1]]))
})

test_that("theta-only E-step preserves exact document modes", {
  case <- make_incremental_case(n=40)
  full <- stm:::estep(
    case$documents, rep(1, 40), FALSE, case$beta, case$lambda,
    case$mu, case$sigma, FALSE
  )
  modes <- stm:::estep.lambda(
    case$documents, rep(1, 40), FALSE, case$beta, case$lambda,
    case$mu, case$sigma, FALSE
  )
  expect_equal(modes, full$lambda, tolerance=1e-12)

  with.bound <- stm:::estep.lambda(
    case$documents, rep(1, 40), FALSE, case$beta, case$lambda,
    case$mu, case$sigma, FALSE, return.bound=TRUE
  )
  expect_equal(with.bound$lambda, full$lambda, tolerance=1e-12)
  expect_equal(with.bound$bound, full$bound, tolerance=1e-12)
})

test_that("memoized M-step equals the batch CTM M-step", {
  case <- make_incremental_case(n=50)
  ss <- stm:::estep(
    case$documents, rep(1, 50), FALSE, case$beta, case$lambda,
    case$mu, case$sigma, FALSE
  )
  settings <- list(
    gamma=list(mode="CTM", maxits=1000),
    sigma=list(prior=.15)
  )
  mapped <- stm:::incremental.mstep(
    ss$beta[[1]], ss$sigma, colSums(ss$lambda), crossprod(ss$lambda),
    NULL, NULL, 50, settings
  )
  expected.mu <- stm:::opt.mu(ss$lambda, "CTM")
  expected.sigma <- stm:::opt.sigma(
    ss$sigma, ss$lambda, expected.mu$mu, settings$sigma$prior
  )
  expected.beta <- stm:::opt.beta(ss$beta, NULL, settings)$beta

  expect_equal(mapped$mu, expected.mu, tolerance=1e-12)
  expect_equal(mapped$sigma, expected.sigma, tolerance=1e-11)
  expect_equal(mapped$beta, expected.beta, tolerance=1e-12)
})

test_that("cached pooled prevalence regression equals the full-data solver", {
  set.seed(29)
  X <- cbind(1, matrix(rnorm(360), 90, 4))
  lambda <- matrix(rnorm(540), 90, 6)
  expected <- stm:::opt.mu(
    lambda, mode="Pooled", covar=X, maxits=1000
  )$gamma
  x.cross <- crossprod(X)
  x.lambda <- crossprod(X, lambda)
  lambda.cross <- crossprod(lambda)
  cached <- sapply(seq_len(ncol(lambda)), function(k) {
    stm:::vb.variational.reg.stats(
      x.lambda[,k], lambda.cross[k,k], x.cross, nrow(X), maxits=1000
    )
  })
  expect_equal(cached, expected, tolerance=1e-12)
})

test_that("incremental defaults adapt block and memory limits", {
  ordinary <- stm:::incremental.defaults(1000, 20, 2632)
  expect_equal(ordinary$batch.size, 160)
  expect_equal(ordinary$nblocks, 7)
  expect_true(ordinary$final.theta)

  large <- stm:::incremental.defaults(
    100000, 50, 50000, max.cache.mb=64
  )
  expect_lte(large$estimated.cache.mb, 64)
  expect_false(large$final.theta)
  expect_equal(large$eval.size, 1000)
  expect_equal(large$early.stop, "both")
  expect_gt(large$batch.size, ordinary$batch.size)
})

test_that("incremental STM is reproducible and theta can be deferred", {
  case <- make_incremental_case(n=60, K=5)
  args <- list(
    documents=case$documents, vocab=case$vocab, K=5,
    init.type="Random", seed=701, max.em.its=3, emtol=0,
    inference="Incremental", verbose=FALSE,
    control=list(incremental.batch.size=20,
                 incremental.final.theta=FALSE)
  )
  first <- do.call(stm, args)
  second <- do.call(stm, args)

  expect_identical(first$beta$logbeta, second$beta$logbeta)
  expect_identical(first$sigma, second$sigma)
  expect_null(first$theta)
  expect_null(first$eta)
  expect_false(first$convergence$theta.computed)

  beta.before <- first$beta$logbeta
  refreshed <- updateTheta(first, case$documents, verbose=FALSE)
  expect_identical(refreshed$beta$logbeta, beta.before)
  expect_equal(dim(refreshed$theta), c(60, 5))
  expect_equal(rowSums(refreshed$theta), rep(1, 60), tolerance=1e-12)
  expect_true(refreshed$convergence$theta.computed)
  expect_true(refreshed$convergence$theta.updated)
  expect_true(is.finite(refreshed$convergence$final.bound))
  expect_true(is.finite(refreshed$convergence$final.bound.per.token))
  expect_gte(refreshed$convergence$theta.update.time, 0)
})

test_that("theta-dependent functions explain how to materialize theta", {
  case <- make_incremental_case(n=30, K=5)
  model <- stm(
    case$documents, case$vocab, K=5, init.type="Random", seed=703,
    max.em.its=1, inference="Incremental", verbose=FALSE,
    control=list(incremental.batch.size=10,
                 incremental.final.theta=FALSE)
  )
  expect_error(findThoughts(model), "updateTheta")
  expect_error(topicCorr(model, verbose=FALSE), "updateTheta")
  expect_equal(
    eval(formals(topicQuality)$labels, list(model=model)),
    seq_len(5)
  )
})

test_that("fresh monitoring bounds are recorded", {
  case <- make_incremental_case(n=60, K=5)
  model <- stm(
    case$documents, case$vocab, K=5, init.type="Random", seed=704,
    max.em.its=3, emtol=0, inference="Incremental", verbose=FALSE,
    control=list(
      incremental.batch.size=20,
      incremental.eval.size=12,
      incremental.eval.every=1,
      incremental.early.stop="both",
      incremental.full.eval.every=1,
      incremental.final.theta=FALSE
    )
  )

  expect_equal(model$convergence$updates, 2 * 3)
  expect_equal(model$convergence$stop.reason, "maximum epochs reached")
  expect_equal(model$convergence$epochs, 3)
  expect_length(model$convergence$eval.bound, 3)
  expect_length(model$convergence$full.bound, 3)
  expect_equal(length(model$settings$incremental$eval.indices), 12)
  expect_true(all(is.finite(model$convergence$eval.per.token)))
  expect_true(all(is.finite(model$convergence$full.per.token)))
  expect_error(
    stm(case$documents, case$vocab, K=5, inference="Incremental",
        verbose=FALSE, control=list(incremental.max.updates=2)),
    "not matched"
  )
})

test_that("incremental restarts use cumulative limits and retain history", {
  case <- make_incremental_case(n=60, K=5)
  common <- list(
    documents=case$documents, vocab=case$vocab, K=5,
    init.type="Random", seed=705, emtol=0,
    inference="Incremental", verbose=FALSE,
    control=list(incremental.batch.size=20,
                 incremental.final.theta=FALSE)
  )
  first <- do.call(stm, c(common, list(max.em.its=2)))
  restarted <- do.call(stm, c(common, list(max.em.its=4, model=first)))

  expect_equal(first$convergence$its, 2)
  expect_equal(restarted$convergence$its, 4)
  expect_equal(restarted$convergence$updates,
               first$convergence$updates + 3)
  expect_equal(head(restarted$convergence$bound, 2),
               first$convergence$bound)
  expect_error(
    do.call(stm, c(common, list(max.em.its=2, model=first))),
    "total epoch limit"
  )
})

test_that("incremental STM supports the default pooled prevalence model", {
  case <- make_incremental_case(n=60, K=5)
  model <- stm(
    case$documents, case$vocab, K=5,
    prevalence=~rating + blog, data=case$meta,
    init.type="Random", seed=702, max.em.its=2, emtol=0,
    inference="Incremental", verbose=FALSE,
    control=list(incremental.batch.size=20,
                 incremental.final.theta=TRUE)
  )
  expect_s3_class(model, "STM")
  expect_equal(nrow(model$theta), 60)
  expect_equal(nrow(model$mu$mu), 4)
  expect_equal(ncol(model$mu$mu), 60)
  expect_false(is.null(model$mu$gamma))
})

test_that("unsupported incremental model components fail clearly", {
  case <- make_incremental_case(n=30, K=5)
  expect_error(
    stm(case$documents, case$vocab, K=5, content=~blog,
        data=case$meta, inference="Incremental", verbose=FALSE),
    "without content covariates"
  )
  expect_error(
    stm(case$documents, case$vocab, K=5,
        prevalence=~rating + blog, data=case$meta,
        gamma.prior="L1", inference="Incremental", verbose=FALSE),
    "default pooled prevalence prior"
  )
  expect_error(
    stm(case$documents, case$vocab, K=5, ngroups=2,
        inference="Incremental", verbose=FALSE),
    "not both"
  )
})
