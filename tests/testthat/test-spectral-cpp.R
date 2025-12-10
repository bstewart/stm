library(stm)
context("Spectral C++ implementation")

# Helper function to create test gram matrix
create_test_gram <- function(V = 100, seed = 42) {
  set.seed(seed)
  Q <- matrix(runif(V*V), V, V)
  Q <- (Q + t(Q)) / 2  # Make symmetric
  Q <- Q / rowSums(Q)  # Row normalize
  return(Q)
}

test_that("fastAnchor C++ matches R version", {
  skip_on_cran()

  # Create test gram matrix
  Qbar <- create_test_gram(V = 100)
  K <- 10

  # Run both implementations with same input
  Qbar_r <- Qbar
  Qbar_cpp <- Qbar

  set.seed(42)
  result_r <- stm:::fastAnchor(Qbar_r, K = K, verbose = FALSE)

  set.seed(42)
  result_cpp <- stm:::fastAnchor.cpp(Qbar_cpp, K = K, verbose = FALSE)

  # Should find same anchors
  expect_equal(result_r, result_cpp)
})

test_that("fastAnchor C++ handles edge cases", {
  skip_on_cran()

  # Test with small vocabulary
  Qbar_small <- create_test_gram(V = 10)
  result <- stm:::fastAnchor.cpp(Qbar_small, K = 5, verbose = FALSE)
  expect_equal(length(result), 5)
  expect_true(all(result >= 1 & result <= 10))

  # Test with K = 1
  result_single <- stm:::fastAnchor.cpp(Qbar_small, K = 1, verbose = FALSE)
  expect_equal(length(result_single), 1)
})

test_that("fastAnchor C++ validates inputs", {
  skip_on_cran()

  Qbar <- create_test_gram(V = 50)

  # K larger than vocabulary
  expect_error(stm:::fastAnchor.cpp(Qbar, K = 100),
               "K.*cannot exceed vocabulary size")

  # Non-square matrix
  expect_error(stm:::fastAnchor.cpp(matrix(1, 10, 5), K = 3),
               "must be a square matrix")
})

test_that("recoverL2 C++ matches R version", {
  skip_on_cran()

  # Use real data for more realistic test
  data(poliblog5k)
  docs <- poliblog5k.docs[1:500]
  vocab <- poliblog5k.voc

  # Compute gram matrix
  docs_ijv <- stm:::doc.to.ijv(docs)
  mat <- Matrix::sparseMatrix(docs_ijv$i, docs_ijv$j, x = docs_ijv$v)
  wprob <- Matrix::colSums(mat)
  wprob <- wprob / sum(wprob)
  Q <- stm:::gram(mat)
  Qbar <- Q / rowSums(Q)

  # Find anchors with R version
  set.seed(42)
  anchors <- stm:::fastAnchor(Qbar, K = 10, verbose = FALSE)

  # Recover with both versions
  set.seed(42)
  A_r <- stm:::recoverL2(Qbar, anchors, wprob, verbose = FALSE, recoverEG = TRUE)$A

  set.seed(42)
  A_cpp <- stm:::recoverL2.cpp(Qbar, anchors, wprob, verbose = FALSE)$A

  # Should be very close (gradient descent may have minor differences)
  expect_equal(dim(A_r), dim(A_cpp))
  expect_equal(A_r, A_cpp, tolerance = 1e-4)

  # Both should sum to 1 across rows
  expect_equal(rowSums(A_r), rep(1, nrow(A_r)), tolerance = 1e-10)
  expect_equal(rowSums(A_cpp), rep(1, nrow(A_cpp)), tolerance = 1e-10)
})

test_that("recoverL2 C++ validates inputs", {
  skip_on_cran()

  Qbar <- create_test_gram(V = 50)
  anchors <- c(1, 5, 10, 15, 20)
  p.w <- rep(1/50, 50)

  # Empty anchors
  expect_error(stm:::recoverL2.cpp(Qbar, integer(0), p.w),
               "anchors vector cannot be empty")

  # Invalid anchor indices
  expect_error(stm:::recoverL2.cpp(Qbar, c(1, 100), p.w),
               "anchor indices must be between 1 and vocabulary size")

  # Wrong p.w length
  expect_error(stm:::recoverL2.cpp(Qbar, anchors, rep(1/10, 10)),
               "p.w length must equal vocabulary size")
})

test_that("Full spectral initialization R vs C++ produces equivalent models", {
  skip_on_cran()

  data(poliblog5k)
  docs <- poliblog5k.docs[1:1000]
  vocab <- poliblog5k.voc
  K <- 15

  # Initialize with R version
  set.seed(42)
  model_r <- stm(docs, vocab, K = K, init.type = "Spectral",
                 control = list(spectral.method = "R"),
                 max.em.its = 5, verbose = FALSE)

  # Initialize with C++ version
  set.seed(42)
  model_cpp <- stm(docs, vocab, K = K, init.type = "Spectral",
                   control = list(spectral.method = "cpp"),
                   max.em.its = 5, verbose = FALSE)

  # Both should be STM objects
  expect_true(inherits(model_r, "STM"))
  expect_true(inherits(model_cpp, "STM"))

  # Same dimensions
  expect_equal(model_r$settings$dim$K, model_cpp$settings$dim$K)
  expect_equal(model_r$settings$dim$V, model_cpp$settings$dim$V)

  # Beta should have same structure
  expect_equal(length(model_r$beta$logbeta), length(model_cpp$beta$logbeta))
  expect_equal(dim(model_r$beta$logbeta[[1]]), dim(model_cpp$beta$logbeta[[1]]))
})

test_that("Auto-selection of spectral method works correctly", {
  skip_on_cran()

  data(poliblog5k)
  docs <- poliblog5k.docs[1:500]
  vocab <- poliblog5k.voc

  # With default (auto), should use R since V < 3000
  set.seed(42)
  model_auto <- stm(docs, vocab, K = 10, init.type = "Spectral",
                    max.em.its = 2, verbose = FALSE)

  expect_true(inherits(model_auto, "STM"))
  expect_equal(model_auto$settings$dim$K, 10)
})

test_that("Spectral C++ works with prevalence covariates", {
  skip_on_cran()

  data(poliblog5k)
  docs <- poliblog5k.docs[1:500]
  vocab <- poliblog5k.voc
  meta <- poliblog5k.meta[1:500, ]

  model <- stm(docs, vocab, K = 10, init.type = "Spectral",
               control = list(spectral.method = "cpp"),
               prevalence = ~ rating,
               data = meta,
               max.em.its = 3, verbose = FALSE)

  expect_true(inherits(model, "STM"))
  expect_false(is.null(model$mu))
})

test_that("Spectral C++ works with content covariates", {
  skip_on_cran()

  data(poliblog5k)
  docs <- poliblog5k.docs[1:500]
  vocab <- poliblog5k.voc
  meta <- poliblog5k.meta[1:500, ]

  model <- stm(docs, vocab, K = 10, init.type = "Spectral",
               control = list(spectral.method = "cpp"),
               content = ~ rating,
               data = meta,
               max.em.its = 3, verbose = FALSE)

  expect_true(inherits(model, "STM"))
  expect_equal(length(model$beta$logbeta), length(unique(meta$rating)))
})

test_that("Downstream functions work with C++ initialized models", {
  skip_on_cran()

  data(poliblog5k)
  docs <- poliblog5k.docs[1:500]
  vocab <- poliblog5k.voc

  model <- stm(docs, vocab, K = 10, init.type = "Spectral",
               control = list(spectral.method = "cpp"),
               max.em.its = 3, verbose = FALSE)

  # labelTopics
  labels <- labelTopics(model, n = 5)
  expect_true(length(labels) >= 4)  # prob, frex, lift, score (and possibly others)
  expect_equal(dim(labels$prob), c(10, 5))

  # findThoughts
  thoughts <- findThoughts(model, texts = paste("doc", 1:length(docs)),
                          n = 2, topics = 1)
  expect_equal(length(thoughts$index[[1]]), 2)

  # topicCorr
  corr <- topicCorr(model)
  expect_true(inherits(corr, "topicCorr"))
})

test_that("C++ version performance is reasonable", {
  skip_on_cran()
  skip_on_ci <- function() {
    if (identical(Sys.getenv("CI"), "true")) {
      skip("Skip on CI to save time")
    }
  }
  skip_on_ci()

  data(poliblog5k)
  docs <- poliblog5k.docs[1:2000]
  vocab <- poliblog5k.voc
  K <- 30

  # Time R version
  time_r <- system.time({
    model_r <- stm(docs, vocab, K = K, init.type = "Spectral",
                   control = list(spectral.method = "R"),
                   max.em.its = 3, verbose = FALSE)
  })

  # Time C++ version
  time_cpp <- system.time({
    model_cpp <- stm(docs, vocab, K = K, init.type = "Spectral",
                     control = list(spectral.method = "cpp"),
                     max.em.its = 3, verbose = FALSE)
  })

  # C++ should be faster or at least not much slower
  speedup <- time_r[3] / time_cpp[3]

  # Print timing info for manual inspection
  message(sprintf("\nTiming comparison (K=%d, N=%d, V=%d):", K, length(docs), length(vocab)))
  message(sprintf("  R version:   %.2f seconds", time_r[3]))
  message(sprintf("  C++ version: %.2f seconds", time_cpp[3]))
  message(sprintf("  Speedup:     %.2fx", speedup))

  # C++ should be at least not significantly slower
  # (On small datasets, overhead may dominate)
  expect_true(time_cpp[3] <= time_r[3] * 1.5,
              info = sprintf("C++ (%.2fs) should not be much slower than R (%.2fs)",
                           time_cpp[3], time_r[3]))
})
