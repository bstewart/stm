make_document_problem <- function(K = 20, V = 400, U = 80, seed = 41) {
  set.seed(seed)
  beta <- matrix(rgamma(K * V, shape = 0.15), nrow = K)
  beta <- beta / rowSums(beta)
  words <- sort(sample.int(V, U))
  counts <- sample.int(8, U, replace = TRUE)
  sigma <- 0.3 ^ abs(outer(seq_len(K - 1), seq_len(K - 1), `-`)) +
    diag(0.2, K - 1)
  sigma_chol <- chol(sigma)
  list(
    eta = rnorm(K - 1, sd = 0.2),
    mu = rnorm(K - 1, sd = 0.1),
    beta = beta[, words, drop = FALSE],
    counts = counts,
    siginv = chol2inv(sigma_chol),
    sigmaentropy = sum(log(diag(sigma_chol)))
  )
}

test_that("zero-copy document kernels preserve numerical results", {
  problem <- make_document_problem(K = 35, U = 120)

  legacy_gradient <- stm:::gradcpp(
    problem$eta, problem$beta, problem$counts,
    problem$mu, problem$siginv
  )
  optimized_gradient <- stm:::gradcpp_nocopy(
    problem$eta, problem$beta, problem$counts,
    problem$mu, problem$siginv
  )
  expect_equal(optimized_gradient, legacy_gradient, tolerance = 1e-10)

  legacy_posterior <- stm:::hpbcpp(
    problem$eta, problem$beta, problem$counts,
    problem$mu, problem$siginv, problem$sigmaentropy
  )
  optimized_posterior <- stm:::hpbcpp_nocopy(
    problem$eta, problem$beta, problem$counts,
    problem$mu, problem$siginv, problem$sigmaentropy
  )
  expect_equal(optimized_posterior$phis, legacy_posterior$phis,
               tolerance = 1e-12)
  expect_equal(optimized_posterior$eta$nu, legacy_posterior$eta$nu,
               tolerance = 1e-10)
  expect_equal(optimized_posterior$bound, legacy_posterior$bound,
               tolerance = 1e-10)
})

test_that("compiled BFGS matches optim at the document mode and bound", {
  problems <- lapply(seq_len(12), function(i) {
    make_document_problem(K = 25, V = 600, U = 40 + 5 * i,
                          seed = 100 + i)
  })
  legacy_bounds <- optimized_bounds <- numeric(length(problems))

  for (i in seq_along(problems)) {
    problem <- problems[[i]]
    legacy <- optim(
      par = problem$eta,
      fn = stm:::lhoodcpp,
      gr = stm:::gradcpp,
      method = "BFGS",
      control = list(maxit = 500),
      beta = problem$beta,
      doc_ct = problem$counts,
      mu = problem$mu,
      siginv = problem$siginv
    )
    optimized <- stm:::optimizecpp(
      eta = problem$eta,
      beta = problem$beta,
      doc_ct = problem$counts,
      mu = problem$mu,
      siginv = problem$siginv,
      control = list(maxit = 500)
    )

    expect_equal(optimized$value, legacy$value, tolerance = 1e-9)
    expect_equal(optimized$par, legacy$par, tolerance = 1e-7)
    expect_equal(optimized$convergence, legacy$convergence)

    legacy_bounds[i] <- stm:::hpbcpp(
      legacy$par, problem$beta, problem$counts,
      problem$mu, problem$siginv, problem$sigmaentropy
    )$bound
    optimized_bounds[i] <- stm:::hpbcpp_nocopy(
      optimized$par, problem$beta, problem$counts,
      problem$mu, problem$siginv, problem$sigmaentropy
    )$bound
  }

  expect_equal(sum(optimized_bounds), sum(legacy_bounds), tolerance = 1e-8)
})

test_that("non-default optimizer features retain the optim fallback", {
  problem <- make_document_problem(K = 10, V = 100, U = 30)
  document <- rbind(seq_len(ncol(problem$beta)), problem$counts)

  fallback <- stm:::logisticnormalcpp(
    eta = problem$eta,
    mu = problem$mu,
    siginv = problem$siginv,
    beta = problem$beta,
    doc = document,
    sigmaentropy = problem$sigmaentropy,
    control = list(maxit = 500, trace = 0)
  )
  expected <- optim(
    par = problem$eta,
    fn = stm:::lhoodcpp,
    gr = stm:::gradcpp,
    method = "BFGS",
    control = list(maxit = 500, trace = 0),
    beta = problem$beta,
    doc_ct = problem$counts,
    mu = problem$mu,
    siginv = problem$siginv
  )

  expect_equal(c(fallback$eta$lambda), expected$par, tolerance = 1e-12)
})
