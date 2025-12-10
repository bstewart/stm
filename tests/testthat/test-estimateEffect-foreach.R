context("estimateEffect() Formula Evaluation in foreach Loops")

library(stm)

# Helper to load the pre-fitted test model
get_test_model <- function() {
  data(gadarianFit, package="stm")
  data(gadarian, package="stm")

  list(model=gadarianFit, metadata=gadarian)
}

test_that("estimateEffect works with variable topics outside foreach", {
  skip_on_cran()

  test_data <- get_test_model()
  model <- test_data$model
  metadata <- test_data$metadata

  # Use variable for topics (the original use case)
  topics <- 1:3
  result <- estimateEffect(topics ~ treatment + s(pid_rep),
                          model, metadata)

  expect_true(inherits(result, "estimateEffect"))
  expect_equal(result$topics, 1:3)
  expect_equal(length(result$parameters), 3)
})

test_that("estimateEffect works with variable topics inside foreach loop", {
  skip_on_cran()
  skip_if_not_installed("foreach")

  library(foreach)

  test_data <- get_test_model()
  model <- test_data$model
  metadata <- test_data$metadata

  # This is the bug from Issue #285 - should now work
  topics <- 1:3
  results <- foreach(i = 1:2) %do% {
    estimateEffect(topics ~ treatment, model, metadata)
  }

  expect_equal(length(results), 2)
  expect_true(all(sapply(results, function(x) inherits(x, "estimateEffect"))))
  expect_true(all(sapply(results, function(x) identical(x$topics, 1:3))))
})

test_that("estimateEffect works in nested function environments", {
  skip_on_cran()
  skip_if_not_installed("foreach")

  library(foreach)

  test_data <- get_test_model()
  model <- test_data$model
  metadata <- test_data$metadata

  # Test nested environments
  test_wrapper <- function() {
    topics <- c(1, 3)  # Different topic selection
    foreach(i = 1:2) %do% {
      estimateEffect(topics ~ treatment, model, metadata)
    }
  }

  results <- test_wrapper()

  expect_equal(length(results), 2)
  expect_true(all(sapply(results, function(x) inherits(x, "estimateEffect"))))
  expect_true(all(sapply(results, function(x) identical(x$topics, c(1, 3)))))
})

test_that("estimateEffect still works with hardcoded topics", {
  skip_on_cran()
  skip_if_not_installed("foreach")

  library(foreach)

  test_data <- get_test_model()
  model <- test_data$model
  metadata <- test_data$metadata

  # Hardcoded topics (the workaround from Issue #285)
  results <- foreach(i = 1:2) %do% {
    estimateEffect(1:3 ~ treatment, model, metadata)
  }

  expect_equal(length(results), 2)
  expect_true(all(sapply(results, function(x) inherits(x, "estimateEffect"))))
  expect_true(all(sapply(results, function(x) identical(x$topics, 1:3))))
})

test_that("estimateEffect gives helpful error for undefined variables", {
  skip_on_cran()

  test_data <- get_test_model()
  model <- test_data$model
  metadata <- test_data$metadata

  # Try to use a variable that doesn't exist
  expect_error(
    estimateEffect(nonexistent_variable ~ treatment, model, metadata),
    "Could not evaluate topic specification"
  )

  # Error message should mention foreach and workarounds
  expect_error(
    estimateEffect(nonexistent_variable ~ treatment, model, metadata),
    "foreach"
  )

  expect_error(
    estimateEffect(nonexistent_variable ~ treatment, model, metadata),
    "Hardcode topic numbers"
  )
})

test_that("estimateEffect works with formula created via as.formula", {
  skip_on_cran()

  test_data <- get_test_model()
  model <- test_data$model
  metadata <- test_data$metadata

  # Test with formula created using as.formula (may not have environment)
  formula_str <- "1:3 ~ treatment"
  formula_obj <- as.formula(formula_str)

  result <- estimateEffect(formula_obj, model, metadata)

  expect_true(inherits(result, "estimateEffect"))
  expect_equal(result$topics, 1:3)
})
