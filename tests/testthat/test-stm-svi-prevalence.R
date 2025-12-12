context("SVI with Prevalence Covariates")

test_that("stm_svi works with simple prevalence formula", {
  skip_on_cran()

  data(poliblog5k, package="stm")
  out <- prepDocuments(poliblog5k.docs, poliblog5k.voc, poliblog5k.meta)

  # Simple formula
  model <- stm_svi(
    documents = out$documents,
    vocab = out$vocab,
    data = out$meta,
    prevalence = ~ rating,
    K = 3,
    max_epochs = 2,
    verbose = FALSE
  )

  expect_true(inherits(model, "STM"))
  expect_false(is.null(model$mu))
  expect_equal(model$settings$dim$K, 3)
  # Should have prevalence covariates
  expect_false(is.null(model$settings$covariates))
})

test_that("stm_svi works with spline prevalence formula", {
  skip_on_cran()

  data(poliblog5k, package="stm")
  out <- prepDocuments(poliblog5k.docs, poliblog5k.voc, poliblog5k.meta)

  # Formula with spline (s() function from mgcv)
  model <- stm_svi(
    documents = out$documents,
    vocab = out$vocab,
    data = out$meta,
    prevalence = ~ s(day),
    K = 3,
    max_epochs = 2,
    verbose = FALSE
  )

  expect_true(inherits(model, "STM"))
  expect_false(is.null(model$mu))
  expect_false(is.null(model$settings$covariates))
})

test_that("stm_svi works with interaction prevalence formula", {
  skip_on_cran()

  data(poliblog5k, package="stm")
  out <- prepDocuments(poliblog5k.docs, poliblog5k.voc, poliblog5k.meta)

  # Formula with interactions
  model <- stm_svi(
    documents = out$documents,
    vocab = out$vocab,
    data = out$meta,
    prevalence = ~ rating * blog,
    K = 3,
    max_epochs = 2,
    verbose = FALSE
  )

  expect_true(inherits(model, "STM"))
  expect_false(is.null(model$mu))
  expect_false(is.null(model$settings$covariates))
})

test_that("stm_svi works with complex prevalence formula", {
  skip_on_cran()

  data(poliblog5k, package="stm")
  out <- prepDocuments(poliblog5k.docs, poliblog5k.voc, poliblog5k.meta)

  # User's original failing case
  model <- stm_svi(
    documents = out$documents,
    vocab = out$vocab,
    data = out$meta,
    prevalence = ~ s(day) + blog + rating,
    K = 5,
    max_epochs = 2,
    verbose = FALSE
  )

  expect_true(inherits(model, "STM"))
  expect_false(is.null(model$mu))
  expect_equal(model$settings$dim$K, 5)
  expect_false(is.null(model$settings$covariates))
})

test_that("stm_svi prevalence matches stm prevalence structure", {
  skip_on_cran()

  data(poliblog5k, package="stm")
  out <- prepDocuments(poliblog5k.docs, poliblog5k.voc, poliblog5k.meta)

  # Same prevalence formula for both
  prev_formula <- ~ rating + blog

  model_stm <- stm(
    documents = out$documents,
    vocab = out$vocab,
    data = out$meta,
    prevalence = prev_formula,
    K = 3,
    max.em.its = 5,
    verbose = FALSE
  )

  model_svi <- stm_svi(
    documents = out$documents,
    vocab = out$vocab,
    data = out$meta,
    prevalence = prev_formula,
    K = 3,
    max_epochs = 2,
    verbose = FALSE
  )

  # Both should have prevalence covariates
  expect_false(is.null(model_stm$settings$covariates))
  expect_false(is.null(model_svi$settings$covariates))

  # Covariate dimensions should match
  expect_equal(
    ncol(model_stm$settings$covariates$X),
    ncol(model_svi$settings$covariates$X)
  )
})

test_that("stm_svi handles intercept-only prevalence correctly", {
  skip_on_cran()

  data(poliblog5k, package="stm")
  out <- prepDocuments(poliblog5k.docs, poliblog5k.voc, poliblog5k.meta)

  # Intercept only (should revert to CTM mode)
  model <- stm_svi(
    documents = out$documents,
    vocab = out$vocab,
    data = out$meta,
    prevalence = ~ 1,
    K = 3,
    max_epochs = 2,
    verbose = FALSE
  )

  expect_true(inherits(model, "STM"))
  # Should have no covariates (intercept-only is CTM)
  expect_null(model$settings$covariates)
})
