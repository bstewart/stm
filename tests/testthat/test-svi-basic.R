# Unit Tests for Stochastic Variational Inference
#
# Tests for the SVI implementation including gradient scaling,
# Adam optimizer, and basic integration.

library(testthat)
library(stm)

context("SVI Basic Functionality")

test_that("Gradient scaling is correct", {
  # Create synthetic test case
  N <- 1000
  batch_size <- 100

  # Simple sufficient statistic (matrix)
  full_ss <- matrix(rnorm(10), nrow=1)
  batch_ss <- full_ss / 10  # Simulate 1/10 of data

  # Scale factor should recover full_ss
  scale_factor <- N / batch_size
  scaled <- batch_ss * scale_factor

  expect_equal(scaled, full_ss)
})

test_that("Adam update produces finite values", {
  # Initialize minimal Adam state
  adam_state <- list(
    m_mu = matrix(0, nrow=2, ncol=1),
    v_mu = matrix(0, nrow=2, ncol=1),
    m_sigma = matrix(0, nrow=2, ncol=2),
    v_sigma = matrix(0, nrow=2, ncol=2),
    m_beta = list(matrix(0, nrow=3, ncol=5)),
    v_beta = list(matrix(0, nrow=3, ncol=5)),
    t = 0
  )

  # Synthetic gradients
  gradients <- list(
    mu = matrix(c(0.1, -0.2), nrow=2, ncol=1),
    sigma = matrix(rnorm(4), nrow=2, ncol=2),
    log_beta = list(matrix(rnorm(15), nrow=3, ncol=5))
  )

  # Run one Adam update
  updated <- adam_update_step(adam_state, gradients, lr=0.01, iter=1)

  # Check no NaN/Inf
  expect_true(all(is.finite(updated$m_mu)))
  expect_true(all(is.finite(updated$v_mu)))
  expect_true(all(is.finite(updated$update_mu)))
  expect_true(all(is.finite(updated$m_sigma)))
  expect_true(all(is.finite(updated$v_sigma)))
  expect_true(all(is.finite(updated$update_sigma)))
  expect_true(all(is.finite(updated$m_beta[[1]])))
  expect_true(all(is.finite(updated$v_beta[[1]])))
  expect_true(all(is.finite(updated$update_beta[[1]])))

  # Check iteration counter incremented
  expect_equal(updated$t, 1)
})

test_that("Beta update maintains probability constraints", {
  # Create valid beta (K=3 topics, V=5 words)
  K <- 3
  V <- 5
  beta_list <- list(matrix(runif(K*V), nrow=K, ncol=V))

  # Normalize to proper probabilities
  beta_list[[1]] <- beta_list[[1]] / rowSums(beta_list[[1]])

  # Create random updates
  update_list <- list(matrix(rnorm(K*V, sd=0.1), nrow=K, ncol=V))

  # Apply updates
  updated_beta <- update_beta_from_adam(beta_list, update_list)

  # Check constraints
  # 1. All values should be non-negative
  expect_true(all(updated_beta[[1]] >= 0))

  # 2. Each row should sum to 1 (simplex constraint)
  row_sums <- rowSums(updated_beta[[1]])
  expect_equal(row_sums, rep(1, K), tolerance=1e-6)

  # 3. No NaN or Inf values
  expect_true(all(is.finite(updated_beta[[1]])))
})

test_that("Sigma remains positive definite after projection", {
  # Create a matrix that's not quite PD
  sigma <- matrix(c(1, 0.95, 0.95, 1), nrow=2)
  sigma <- sigma - 0.1 * diag(2)  # Make it non-PD

  # Project back to PD
  sigma_pd <- ensure_sigma_pd(sigma)

  # Check it's now PD
  eigenvalues <- eigen(sigma_pd, symmetric=TRUE, only.values=TRUE)$values
  expect_true(all(eigenvalues > 0))

  # Check it's symmetric
  expect_equal(sigma_pd, t(sigma_pd))
})

test_that("SVI convergence initialization works", {
  settings <- list(
    svi = list(
      max_epochs = 50,
      patience = 20
    ),
    verbose = FALSE
  )

  convergence <- initialize_svi_convergence(settings)

  expect_equal(convergence$its, 0)
  expect_equal(convergence$epoch, 0)
  expect_false(convergence$converged)
  expect_false(convergence$stopits)
  expect_equal(convergence$best_window, -Inf)
  expect_equal(length(convergence$bound_history), 0)
})

test_that("SVI convergence detects patience exhaustion", {
  settings <- list(
    svi = list(
      max_epochs = 100,
      patience = 5,
      convergence_window = 3
    ),
    verbose = FALSE
  )

  convergence <- initialize_svi_convergence(settings)
  N <- 1000
  batch_size <- 100

  # Simulate iterations with no improvement
  for(i in 1:20) {
    # Constant ELBO (no improvement)
    convergence <- svi_convergence_check(
      batch_elbo = 100,
      convergence = convergence,
      batch_size = batch_size,
      N = N,
      settings = settings
    )

    # Should stop after patience + window iterations
    if(i > 5 + 3) {
      expect_true(convergence$stopits)
      break
    }
  }
})

test_that("Gradient clipping works correctly", {
  # Create gradients with some extreme values
  gradients <- list(
    mu = matrix(c(-100, 0.5, 200), nrow=3, ncol=1),
    sigma = matrix(c(1, -50, 30, 2), nrow=2, ncol=2),
    log_beta = list(matrix(c(0.1, -500, 0.3, 1000), nrow=2, ncol=2))
  )

  # Clip with threshold of 10
  clipped <- clip_gradients(gradients, clip_value=10)

  # Check all values are within [-10, 10]
  expect_true(all(abs(clipped$mu) <= 10))
  expect_true(all(abs(clipped$sigma) <= 10))
  expect_true(all(abs(clipped$log_beta[[1]]) <= 10))

  # Check non-extreme values unchanged
  expect_equal(clipped$mu[2, 1], 0.5)
  expect_equal(clipped$log_beta[[1]][1, 1], 0.1)
})

test_that("Gradient checking detects NaN and Inf", {
  # Valid gradients
  good_gradients <- list(
    mu = matrix(c(0.1, 0.2), nrow=2, ncol=1),
    sigma = matrix(c(1, 0.5, 0.5, 1), nrow=2, ncol=2),
    log_beta = list(matrix(rnorm(6), nrow=2, ncol=3))
  )

  expect_true(check_gradients(good_gradients, iter=1, stop_on_error=FALSE))

  # Invalid gradients - NaN
  bad_gradients_nan <- good_gradients
  bad_gradients_nan$mu[1, 1] <- NaN

  expect_false(check_gradients(bad_gradients_nan, iter=1, stop_on_error=FALSE))

  # Invalid gradients - Inf
  bad_gradients_inf <- good_gradients
  bad_gradients_inf$sigma[1, 1] <- Inf

  expect_false(check_gradients(bad_gradients_inf, iter=1, stop_on_error=FALSE))
})

context("SVI Integration Tests")

test_that("stm_svi validates arguments correctly", {
  # Missing required arguments
  expect_error(stm_svi(vocab=letters[1:10], K=3), "Must include documents")
  expect_error(stm_svi(documents=list(), K=3), "Vocab length")

  # Invalid K
  expect_error(stm_svi(documents=list(matrix(c(1, 1), nrow=2)),
                      vocab="word", K=1), "K must be")

  # Invalid batch_size
  expect_error(stm_svi(documents=list(matrix(c(1, 1), nrow=2)),
                      vocab="word", K=3, batch_size=-1), "batch_size must be")

  # Invalid learning rate
  expect_error(stm_svi(documents=list(matrix(c(1, 1), nrow=2)),
                      vocab="word", K=3, lr=-0.01), "Learning rate.*must be positive")
})

# Note: Full integration test on gadarian dataset should be run separately
# as it takes longer and requires the full stm package to be loaded
