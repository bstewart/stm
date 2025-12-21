# Test adaptive hyperparameter selection for stm_svi()

test_that("Adaptive defaults for small corpus", {
  skip_on_cran()

  # Use gadarian data (small corpus)
  data(gadarian)
  docs <- gadarian$documents[1:1000]
  vocab <- gadarian$vocab

  # Should work without specifying any hyperparameters
  model <- suppressWarnings(
    stm_svi(docs, vocab, K=20, verbose=FALSE)
  )

  expect_true(inherits(model, "STM"))

  # Should have used adaptive batch_size around 80 (20% of 1000, or 4*K=80)
  # Allow some tolerance
  expect_true(model$settings$svi$batch_size >= 60 &&
              model$settings$svi$batch_size <= 100)

  # Should have set high max_epochs for small corpus
  expect_true(model$settings$svi$max_epochs >= 10)

  # Should have set lr
  expect_true(!is.null(model$settings$svi$lr))
  expect_true(model$settings$svi$lr > 0 && model$settings$svi$lr < 0.1)

  # Should have set patience and convergence_window
  expect_true(!is.null(model$settings$svi$patience))
  expect_true(!is.null(model$settings$svi$convergence_window))
})

test_that("Adaptive defaults for medium corpus with prevalence", {
  skip_on_cran()
  skip_if_not_installed("poliblog5k")

  # Use poliblog5k data
  data(poliblog5k)

  model <- suppressWarnings(
    stm_svi(poliblog5k.docs, poliblog5k.voc, K=50,
           prevalence=~rating, data=poliblog5k.meta,
           verbose=FALSE)
  )

  expect_true(inherits(model, "STM"))

  # Should have set gamma_update_every for prevalence
  expect_true(!is.null(model$settings$svi$gamma_update_every))
  expect_equal(model$settings$svi$gamma_update_every, 1)

  # Batch size should scale with K (4*50 = 200)
  expect_true(model$settings$svi$batch_size >= 150 &&
              model$settings$svi$batch_size <= 250)
})

test_that("User-specified parameters override adaptive defaults", {
  skip_on_cran()

  # Use gadarian data
  data(gadarian)
  docs <- gadarian$documents[1:1000]
  vocab <- gadarian$vocab

  # Specify all parameters explicitly
  model <- suppressWarnings(
    stm_svi(docs, vocab, K=20,
           batch_size=256, lr=0.05, max_epochs=10, patience=5,
           convergence_window=15,
           verbose=FALSE)
  )

  # Should use exact user values, not adaptive
  expect_equal(model$settings$svi$batch_size, min(256, 1000))  # Can't exceed N
  expect_equal(model$settings$svi$lr, 0.05)
  expect_equal(model$settings$svi$max_epochs, 10)
  expect_equal(model$settings$svi$patience, 5)
  expect_equal(model$settings$svi$convergence_window, 15)
})

test_that("Adaptive defaults scale appropriately with K", {
  skip_on_cran()

  # Use gadarian data
  data(gadarian)
  docs <- gadarian$documents[1:2000]
  vocab <- gadarian$vocab

  # Test with K=10
  model_k10 <- suppressWarnings(
    stm_svi(docs, vocab, K=10, verbose=FALSE, max_iters=3)
  )

  # Test with K=50
  model_k50 <- suppressWarnings(
    stm_svi(docs, vocab, K=50, verbose=FALSE, max_iters=3)
  )

  # Batch size should increase with K (4*K base)
  expect_true(model_k50$settings$svi$batch_size > model_k10$settings$svi$batch_size)

  # Learning rate should be reduced for high K (K>75 gets 0.8x multiplier)
  # K=50 shouldn't trigger this, so lr should be similar
  expect_true(abs(model_k10$settings$svi$lr - model_k50$settings$svi$lr) < 0.01)
})

test_that("Batch size validation works correctly", {
  skip_on_cran()

  # Use small corpus
  data(gadarian)
  docs <- gadarian$documents[1:100]
  vocab <- gadarian$vocab

  # Adaptive should not exceed N
  model <- suppressWarnings(
    stm_svi(docs, vocab, K=20, verbose=FALSE, max_iters=3)
  )

  expect_true(model$settings$svi$batch_size <= 100)
})

test_that("Compute_svi_defaults produces valid outputs", {
  # Test the helper function directly
  defaults <- stm:::compute_svi_defaults(N=5000, K=50, V=2000, has_prevalence=FALSE)

  expect_true(is.list(defaults))
  expect_true(all(c("batch_size", "lr", "max_epochs", "patience",
                    "convergence_window", "iters_per_epoch",
                    "target_total_iters", "gradient_scale") %in% names(defaults)))

  # Check ranges
  expect_true(defaults$batch_size >= 32 && defaults$batch_size <= 1024)
  expect_true(defaults$lr >= 0.005 && defaults$lr <= 0.02)
  expect_true(defaults$max_epochs >= 5 && defaults$max_epochs <= 100)
  expect_true(defaults$patience >= 10 && defaults$patience <= 50)
  expect_true(defaults$convergence_window >= 5 && defaults$convergence_window <= 20)

  # With prevalence
  defaults_prev <- stm:::compute_svi_defaults(N=5000, K=50, V=2000, has_prevalence=TRUE)
  expect_true(!is.null(defaults_prev$gamma_update_every))
  expect_equal(defaults_prev$gamma_update_every, 1)

  # Without prevalence
  expect_true(is.null(defaults$gamma_update_every))
})
