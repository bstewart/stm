make_spectral_test_documents <- function() {
  list(
    rbind(as.integer(c(1, 2, 4)), as.integer(c(4, 2, 1))),
    rbind(as.integer(c(1, 3, 5)), as.integer(c(1, 4, 2))),
    rbind(as.integer(c(2, 3, 6)), as.integer(c(3, 1, 3))),
    rbind(as.integer(c(1, 4, 6)), as.integer(c(2, 3, 2))),
    rbind(as.integer(c(2, 5, 6)), as.integer(c(1, 2, 4))),
    rbind(as.integer(c(3, 4, 5)), as.integer(c(3, 2, 1)))
  )
}

test_that("native exact spectral kernels match the R reference", {
  documents <- make_spectral_test_documents()
  vocabulary_size <- 6L
  triplet <- stm:::doc.to.ijv(documents)
  document_term <- Matrix::sparseMatrix(
    i=triplet$i, j=triplet$j, x=triplet$v,
    dims=c(length(documents), vocabulary_size))

  reference_q <- stm:::gram(document_term)
  reference_probability <- Matrix::colSums(document_term)
  reference_probability <- reference_probability/sum(reference_probability)
  native <- stm:::spectralGramCpp(documents, vocabulary_size)

  expect_equal(native$Q, reference_q, tolerance=1e-14)
  expect_equal(native$wprob, as.numeric(reference_probability), tolerance=0)

  qbar <- reference_q/rowSums(reference_q)
  reference_anchor <- stm:::fastAnchor(qbar, K=3, verbose=FALSE)
  native_anchor <- stm:::fastAnchorCpp(qbar, topics=3)
  expect_identical(native_anchor, as.integer(reference_anchor))

  reference_beta <- stm:::recoverL2(
    qbar, reference_anchor, reference_probability,
    verbose=FALSE, recoverEG=TRUE)$A
  native_beta <- stm:::recoverL2Cpp(
    qbar, native_anchor, reference_probability,
    tolerance=1e-7, maximum_iterations=500)$A
  expect_equal(native_beta, reference_beta, tolerance=1e-13)
})

test_that("streamed projection has exact normalizers and reproducible hashes", {
  documents <- make_spectral_test_documents()
  vocabulary_size <- 6L
  exact <- stm:::spectralGramCpp(documents, vocabulary_size)
  set.seed(1)
  first <- stm:::spectralProjectionCpp(
    documents, vocabulary_size, projection_dimension=5,
    projection_sparsity=.4, seed=9)
  set.seed(999)
  second <- stm:::spectralProjectionCpp(
    documents, vocabulary_size, projection_dimension=5,
    projection_sparsity=.4, seed=9)

  expect_equal(first$row_normalizer, rowSums(exact$Q), tolerance=1e-13)
  expect_identical(first$Qbar, second$Qbar)
  expect_identical(first$wprob, second$wprob)
  expect_true(all(is.finite(first$Qbar)))
})

test_that("original-space recovery is exact for a fixed anchor set", {
  documents <- make_spectral_test_documents()
  vocabulary_size <- 6L
  exact <- stm:::spectralGramCpp(documents, vocabulary_size)
  qbar <- exact$Q/rowSums(exact$Q)
  anchors <- stm:::fastAnchorCpp(qbar, topics=3)
  reference <- stm:::recoverL2Cpp(
    qbar, anchors, exact$wprob,
    tolerance=1e-7, maximum_iterations=500)$A
  streamed <- stm:::recoverL2OriginalCpp(
    documents, vocabulary_size, anchors,
    row_normalizer=rowSums(exact$Q),
    word_probability=exact$wprob,
    tolerance=1e-7, maximum_iterations=500)$A

  expect_equal(streamed, reference, tolerance=1e-13)
  expect_equal(rowSums(streamed), rep(1, 3), tolerance=1e-13)
})

test_that("exact recovery takes over when exponentiated gradient does not converge", {
  # Short, peaked documents make the row-normalized moment rows concentrated,
  # so the fixed exponentiated-gradient step of 50 oscillates.  The solver must
  # then return the exact simplex-constrained least-squares solution.
  temp <- textProcessor(gadarian$open.ended.response, metadata=gadarian,
                        verbose=FALSE)
  out <- prepDocuments(temp$documents, temp$vocab, temp$meta, verbose=FALSE)
  V <- length(out$vocab); K <- 3
  moment <- stm:::spectralGramCpp(out$documents, V)
  qbar <- moment$Q/rowSums(moment$Q)
  anchors <- stm:::fastAnchorCpp(qbar, K)
  recovered <- stm:::recoverL2Cpp(qbar, anchors, moment$wprob,
                                  tolerance=1e-7, maximum_iterations=500L)
  expect_true(any(recovered$fallback))
  expect_true(all(recovered$converged))
  expect_equal(rowSums(recovered$A), rep(1, K), tolerance=1e-12)
  expect_equal(rowSums(recovered$weights), rep(1, V), tolerance=1e-12)

  # Fallback words agree with quadprog's exact solution of the same problem,
  # and words where exponentiated gradient converged keep the legacy solution.
  X <- qbar[anchors,,drop=FALSE]; XtX <- tcrossprod(X)
  Amat <- cbind(1, diag(1, K)); bvec <- c(1, rep(0, K))
  words <- head(which(recovered$fallback), 25)
  for(i in words) {
    qp <- quadprog::solve.QP(Dmat=XtX, dvec=as.numeric(X %*% qbar[i,]),
                             Amat=Amat, bvec=bvec, meq=1)$solution
    qp[qp < 0] <- 0
    expect_equal(recovered$weights[i,], qp/sum(qp), tolerance=1e-6)
  }
  legacy.words <- which(!recovered$fallback & !(seq_len(V) %in% anchors))
  for(i in head(legacy.words, 3)) {
    legacy <- stm:::expgrad(X, qbar[i,], XtX)
    expect_true(legacy$converged)
    expect_equal(recovered$weights[i,], legacy$par, tolerance=1e-12)
  }
})

test_that("recoverL2 dispatches exponentiated gradient to compiled code", {
  documents <- make_spectral_test_documents()
  exact <- stm:::spectralGramCpp(documents, 6L)
  qbar <- exact$Q/rowSums(exact$Q)
  anchors <- stm:::fastAnchorCpp(qbar, topics=3)
  via.r <- stm:::recoverL2(qbar, anchors, exact$wprob, verbose=FALSE)$A
  native <- stm:::recoverL2Cpp(qbar, anchors, exact$wprob)$A
  expect_identical(via.r, native)
})

test_that("the restricted moment matches gram() on the truncated matrix", {
  documents <- make_spectral_test_documents()
  vocabulary_size <- 6L
  keep <- c(5L, 2L, 6L, 1L) # most frequent first, as maxV orders them
  triplet <- stm:::doc.to.ijv(documents)
  document_term <- Matrix::sparseMatrix(
    i=triplet$i, j=triplet$j, x=triplet$v,
    dims=c(length(documents), vocabulary_size))
  reference <- stm:::gram(document_term[, keep])
  native <- stm:::spectralGramCpp(documents, vocabulary_size, keep=keep)
  expect_equal(native$Q, reference, tolerance=1e-14)
  expect_equal(length(native$wprob), vocabulary_size)
  expect_equal(native$wprob, as.numeric(Matrix::colSums(document_term)/sum(document_term)))
})

test_that("the sketch is opt-in and the exact path is the default", {
  skip_on_cran()
  set.seed(5)
  docs <- poliblog5k.docs[1:400]
  prepared <- prepDocuments(docs, poliblog5k.voc, lower.thresh=0, verbose=FALSE)
  exact <- stm(prepared$documents, prepared$vocab, K=12, max.em.its=0,
               verbose=FALSE)
  expect_false(exact$settings$init$projected)
  expect_equal(exact$settings$init$p, 1024)
  sketched <- stm(prepared$documents, prepared$vocab, K=12, max.em.its=0,
                  verbose=FALSE, control=list(SpectralRP=TRUE, rp.p=128))
  expect_true(sketched$settings$init$projected)
  expect_equal(dim(sketched$beta$logbeta[[1]]), dim(exact$beta$logbeta[[1]]))
  expect_equal(rowSums(exp(sketched$beta$logbeta[[1]])), rep(1, 12), tolerance=1e-10)
  # the sketch is deterministic and independent of the R seed
  set.seed(99)
  again <- stm(prepared$documents, prepared$vocab, K=12, max.em.its=0,
               verbose=FALSE, control=list(SpectralRP=TRUE, rp.p=128))
  expect_identical(again$beta$logbeta, sketched$beta$logbeta)
  # maxV runs through the native path and reintroduces dropped words
  truncated <- stm(prepared$documents, prepared$vocab, K=12, max.em.its=0,
                   verbose=FALSE, control=list(maxV=300))
  expect_equal(rowSums(exp(truncated$beta$logbeta[[1]])), rep(1, 12), tolerance=1e-10)
  expect_true(all(is.finite(truncated$beta$logbeta[[1]])))
})

test_that("recover.exact solves every word exactly", {
  documents <- make_spectral_test_documents()
  exact <- stm:::spectralGramCpp(documents, 6L)
  qbar <- exact$Q/rowSums(exact$Q)
  anchors <- stm:::fastAnchorCpp(qbar, topics=3)
  all.exact <- stm:::recoverL2Cpp(qbar, anchors, exact$wprob, exact_only=TRUE)
  expect_true(all(all.exact$fallback[-anchors]))
  expect_true(all(all.exact$converged))
  expect_equal(rowSums(all.exact$A), rep(1, 3), tolerance=1e-12)
  hybrid <- stm:::recoverL2Cpp(qbar, anchors, exact$wprob)
  # the exact optimum has squared error no larger than the EG solution
  X <- qbar[anchors,,drop=FALSE]; XtX <- tcrossprod(X)
  sse <- function(w, i) sum((as.numeric(qbar[i,] %*% t(X)) - as.numeric(w %*% XtX))^2)
  for(i in setdiff(1:6, anchors)) {
    expect_lte(sse(all.exact$weights[i,], i), sse(hybrid$weights[i,], i) + 1e-12)
  }
  skip_on_cran()
  prepared <- prepDocuments(poliblog5k.docs[1:200], poliblog5k.voc, lower.thresh=0, verbose=FALSE)
  fit <- stm(prepared$documents, prepared$vocab, K=5, max.em.its=0, verbose=FALSE,
             control=list(recover.exact=TRUE))
  expect_true(fit$settings$init$recover.exact)
  expect_equal(rowSums(exp(fit$beta$logbeta[[1]])), rep(1, 5), tolerance=1e-10)
})
