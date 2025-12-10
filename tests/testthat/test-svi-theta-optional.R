context("SVI Optional Theta")

# Helper to create small test dataset
get_test_data <- function() {
  data(poliblog5k, package="stm")
  # Use subset and prepDocuments to ensure sequential word indices
  docs_subset <- poliblog5k.docs[1:200]
  prep <- prepDocuments(docs_subset, poliblog5k.voc, lower.thresh=0)
  list(docs=prep$documents, vocab=prep$vocab)
}

test_that("stm_svi with compute_final_theta=FALSE returns NULL theta", {
  skip_on_cran()

  test_data <- get_test_data()
  docs <- test_data$docs
  vocab <- test_data$vocab

  # Fit without theta
  model <- stm_svi(docs, vocab, K=3,
                   compute_final_theta=FALSE,
                   max_epochs=2,
                   init.type="Random",
                   seed=123,
                   verbose=FALSE)

  expect_null(model$theta)
  expect_null(model$eta)
  expect_false(model$convergence$theta_computed)
  expect_true(model$svi)
})

test_that("stm_svi with compute_final_theta=TRUE returns theta", {
  skip_on_cran()

  test_data <- get_test_data()
  docs <- test_data$docs
  vocab <- test_data$vocab

  # Fit with theta
  model <- stm_svi(docs, vocab, K=3,
                   compute_final_theta=TRUE,
                   max_epochs=2,
                   init.type="Random",
                   seed=123,
                   verbose=FALSE)

  expect_false(is.null(model$theta))
  expect_false(is.null(model$eta))
  expect_true(model$convergence$theta_computed)
  expect_equal(nrow(model$theta), length(docs))
  expect_equal(ncol(model$theta), 3)
})

test_that("update_svi_theta computes theta for model without it", {
  skip_on_cran()

  test_data <- get_test_data()
  docs <- test_data$docs
  vocab <- test_data$vocab

  # Fit without theta
  model <- stm_svi(docs, vocab, K=3,
                   compute_final_theta=FALSE,
                   max_epochs=2,
                   init.type="Random",
                   seed=123,
                   verbose=FALSE)

  expect_null(model$theta)

  # Compute theta
  model <- update_svi_theta(model, docs, verbose=FALSE)

  expect_false(is.null(model$theta))
  expect_false(is.null(model$eta))
  expect_true(model$convergence$theta_computed)
  expect_equal(nrow(model$theta), length(docs))
})

test_that("Functions requiring theta raise helpful errors", {
  skip_on_cran()

  test_data <- get_test_data()
  docs <- test_data$docs
  vocab <- test_data$vocab

  # Fit without theta
  model <- stm_svi(docs, vocab, K=3,
                   compute_final_theta=FALSE,
                   max_epochs=2,
                   init.type="Random",
                   seed=123,
                   verbose=FALSE)

  # These should all error with helpful message
  expect_error(findThoughts(model, n=3), "update_svi_theta")
  expect_error(plot(model, type="summary"), "update_svi_theta")
  expect_error(topicCorr(model), "update_svi_theta")

  # After computing theta, should work
  model <- update_svi_theta(model, docs, verbose=FALSE)

  expect_silent(thoughts <- findThoughts(model, n=3))
  expect_silent(plot(model, type="summary"))
  expect_silent(corr <- topicCorr(model))
})

test_that("topicQuality works without theta using settings$dim$K", {
  skip_on_cran()

  test_data <- get_test_data()
  docs <- test_data$docs
  vocab <- test_data$vocab

  # Fit without theta
  model <- stm_svi(docs, vocab, K=3,
                   compute_final_theta=FALSE,
                   max_epochs=2,
                   init.type="Random",
                   seed=123,
                   verbose=FALSE)

  # topicQuality should work because it uses settings$dim$K
  expect_silent(topicQuality(model, docs))
})
