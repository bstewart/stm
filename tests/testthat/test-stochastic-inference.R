make_stochastic_case <- function(n=1000, seed=31) {
  prepared <- prepDocuments(
    poliblog5k.docs[seq_len(n)], poliblog5k.voc,
    poliblog5k.meta[seq_len(n),,drop=FALSE],
    lower.thresh=0, verbose=FALSE
  )
  list(documents=prepared$documents, vocab=prepared$vocab, meta=prepared$meta)
}

test_that("stochastic inference runs its two phases and is reproducible", {
  case <- make_stochastic_case()
  args <- list(
    documents=case$documents, vocab=case$vocab, K=5,
    init.type="Random", seed=811, max.em.its=4, emtol=0,
    inference="Stochastic", verbose=FALSE,
    control=list(stochastic.batch.size=50, stochastic.max.epochs=2,
                 stochastic.tol=0)
  )
  first <- do.call(stm, args)
  second <- do.call(stm, args)
  expect_s3_class(first, "STM")
  expect_identical(first$inference, "Stochastic")
  expect_identical(first$beta$logbeta, second$beta$logbeta)
  expect_identical(first$theta, second$theta)

  record <- first$convergence$stochastic
  expect_equal(record$epochs, 2)
  expect_equal(record$updates, 2 * 20)
  expect_equal(record$batch.size, 50)
  expect_length(record$minibatch.bound, 40)
  expect_true(all(is.finite(record$minibatch.bound)))
  expect_true(all(is.finite(record$checkpoint.bound)))
  expect_equal(record$stop.reason, "maximum stochastic epochs reached")
  # the memoized phase reused the minibatches as blocks and its cache
  expect_equal(first$convergence$blocks, 20)
  expect_equal(first$convergence$batch.size, 50)
  expect_equal(first$convergence$its, 4)
  expect_equal(dim(first$theta), c(1000, 5))
  expect_equal(rowSums(first$theta), rep(1, 1000), tolerance=1e-12)
})

test_that("the stochastic phase supports prevalence covariates and can stop early", {
  case <- make_stochastic_case()
  model <- stm(
    case$documents, case$vocab, K=5, prevalence=~rating + s(day),
    data=case$meta, init.type="Random", seed=812, max.em.its=3, emtol=0,
    inference="Stochastic", verbose=FALSE,
    control=list(stochastic.batch.size=50, stochastic.max.epochs=5,
                 stochastic.tol=1) # any gain counts as flat: stops at the earliest allowed point
  )
  expect_false(is.null(model$mu$gamma))
  expect_equal(ncol(model$mu$mu), 1000)
  expect_equal(model$convergence$stochastic$stop.reason,
               "bound gain per pass fell below stochastic.tol")
  expect_lt(model$convergence$stochastic$epochs, 2)
})

test_that("stochastic.finish=FALSE returns the stochastic solution", {
  case <- make_stochastic_case()
  model <- stm(
    case$documents, case$vocab, K=5, init.type="Random", seed=813,
    inference="Stochastic", verbose=FALSE,
    control=list(stochastic.batch.size=50, stochastic.max.epochs=1,
                 stochastic.tol=0, stochastic.finish=FALSE)
  )
  expect_identical(model$inference, "Stochastic")
  expect_equal(model$convergence$its, 0)
  expect_equal(model$convergence$epochs, 1)
  expect_true(is.finite(model$convergence$final.bound))
  expect_equal(dim(model$theta), c(1000, 5))
  deferred <- stm(
    case$documents, case$vocab, K=5, init.type="Random", seed=813,
    inference="Stochastic", verbose=FALSE,
    control=list(stochastic.batch.size=50, stochastic.max.epochs=1,
                 stochastic.tol=0, stochastic.finish=FALSE,
                 incremental.final.theta=FALSE)
  )
  expect_null(deferred$theta)
  expect_identical(deferred$beta$logbeta, model$beta$logbeta)
  refreshed <- updateTheta(deferred, case$documents, verbose=FALSE)
  # the deferred pass starts each document from the prior mean rather than its
  # last stochastic iterate, so agreement is at the optimizer's tolerance
  expect_equal(refreshed$theta, model$theta, tolerance=1e-3)
  expect_equal(refreshed$convergence$final.bound, model$convergence$final.bound,
               tolerance=1e-6)
})

test_that("small corpora skip the stochastic phase", {
  case <- make_stochastic_case(n=120)
  model <- stm(case$documents, case$vocab, K=4, init.type="Random", seed=814,
               max.em.its=2, emtol=0, inference="Stochastic", verbose=FALSE)
  expect_equal(model$convergence$stochastic$stop.reason, "too few minibatches")
  expect_equal(model$convergence$stochastic$updates, 0)
  expect_identical(model$inference, "Stochastic")
  expect_equal(model$convergence$its, 2)
})

test_that("stochastic inference validates its inputs", {
  case <- make_stochastic_case(n=200)
  expect_error(
    stm(case$documents, case$vocab, K=4, inference="Stochastic",
        verbose=FALSE, control=list(stochastic.kappa=0.4)),
    "kappa"
  )
  expect_error(
    stm(case$documents, case$vocab, K=4, inference="Stochastic",
        verbose=FALSE, control=list(stochastic.finish="yes")),
    "finish"
  )
  fit <- stm(case$documents, case$vocab, K=4, init.type="Random", seed=1,
             max.em.its=1, inference="Stochastic", verbose=FALSE)
  expect_error(
    stm(case$documents, case$vocab, K=4, inference="Stochastic",
        model=fit, verbose=FALSE),
    "Restarting is not supported"
  )
  expect_error(
    stm(case$documents, case$vocab, K=4, content=~rating, data=case$meta,
        inference="Stochastic", verbose=FALSE),
    "without content covariates"
  )
})
