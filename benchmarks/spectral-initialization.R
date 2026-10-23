#!/usr/bin/env Rscript

# Phase-level benchmark for spectral initialization. The four rows are
# cumulative ablations: legacy R, exact native kernels, streamed anchor search
# with strict original-space recovery, and the final recovery tolerance.

suppressPackageStartupMessages(library(stm))

args <- commandArgs(trailingOnly=TRUE)
repetitions <- if(length(args)) as.integer(args[[1]]) else 3L
output <- if(length(args) >= 2) args[[2]] else NULL
if(is.na(repetitions) || repetitions < 1) stop("repetitions must be positive")

data(poliblog5k, package="stm")

prepare_case <- function(document_count, topics) {
  if(document_count == length(poliblog5k.docs)) {
    return(list(documents=poliblog5k.docs, vocab=poliblog5k.voc, K=topics))
  }
  prepared <- prepDocuments(
    poliblog5k.docs[seq_len(document_count)], poliblog5k.voc,
    lower.thresh=0, verbose=FALSE)
  list(documents=prepared$documents, vocab=prepared$vocab, K=topics)
}

cases <- list(
  poliblog_300=prepare_case(300, 10),
  poliblog_1000=prepare_case(1000, 20),
  poliblog_5000=prepare_case(5000, 30)
)

restore_words <- function(beta, valid, vocabulary_size) {
  if(all(valid)) return(beta)
  restored <- matrix(0, nrow=nrow(beta), ncol=vocabulary_size)
  if(ncol(beta) == sum(valid)) restored[,valid] <- beta else restored <- beta
  restored <- restored + .001/vocabulary_size
  restored/rowSums(restored)
}

legacy_exact <- function(case) {
  triplet <- stm:::doc.to.ijv(case$documents)
  document_term <- Matrix::sparseMatrix(
    triplet$i, triplet$j, x=triplet$v,
    dims=c(length(case$documents), length(case$vocab)))
  word_probability <- Matrix::colSums(document_term)
  word_probability <- word_probability/sum(word_probability)

  moment_time <- system.time(moment <- stm:::gram(document_term))[["elapsed"]]
  valid <- rowSums(moment) > 0
  qbar <- moment[valid,valid,drop=FALSE]/rowSums(moment)[valid]
  anchor_time <- system.time(
    anchor <- stm:::fastAnchor(qbar, case$K, verbose=FALSE))[["elapsed"]]
  recovery_time <- system.time(
    beta <- stm:::recoverL2(
      qbar, anchor, word_probability[valid], verbose=FALSE)$A)[["elapsed"]]
  beta <- restore_words(beta, valid, length(case$vocab))
  list(
    moment=moment_time, anchor=anchor_time, recovery=recovery_time,
    total=moment_time + anchor_time + recovery_time,
    anchors=which(valid)[anchor], beta=beta)
}

native_exact <- function(case) {
  moment_time <- system.time(
    moment <- stm:::spectralGramCpp(case$documents, length(case$vocab))
  )[["elapsed"]]
  valid <- rowSums(moment$Q) > 0
  qbar <- moment$Q[valid,valid,drop=FALSE]/rowSums(moment$Q)[valid]
  anchor_time <- system.time(
    anchor <- stm:::fastAnchorCpp(qbar, case$K))[["elapsed"]]
  recovery_time <- system.time(
    beta <- stm:::recoverL2Cpp(
      qbar, anchor, moment$wprob[valid],
      tolerance=1e-7, maximum_iterations=500)$A)[["elapsed"]]
  beta <- restore_words(beta, valid, length(case$vocab))
  list(
    moment=moment_time, anchor=anchor_time, recovery=recovery_time,
    total=moment_time + anchor_time + recovery_time,
    anchors=which(valid)[anchor], beta=beta)
}

projected_original <- function(case, tolerance) {
  dimension <- max(256L, case$K)
  moment_time <- system.time(
    sketch <- stm:::spectralProjectionCpp(
      case$documents, length(case$vocab), dimension,
      projection_sparsity=.05, seed=9L))[["elapsed"]]
  valid <- sketch$row_normalizer > 0
  anchor_time <- system.time(
    anchor_reduced <- stm:::fastAnchorCpp(
      sketch$Qbar[valid,,drop=FALSE], case$K))[["elapsed"]]
  anchors <- which(valid)[anchor_reduced]
  recovery_time <- system.time(
    beta <- stm:::recoverL2OriginalCpp(
      case$documents, length(case$vocab), anchors,
      sketch$row_normalizer, sketch$wprob,
      tolerance=tolerance, maximum_iterations=500)$A)[["elapsed"]]
  beta <- restore_words(beta, valid, length(case$vocab))
  list(
    moment=moment_time, anchor=anchor_time, recovery=recovery_time,
    total=moment_time + anchor_time + recovery_time,
    anchors=anchors, beta=beta)
}

methods <- list(
  legacy=legacy_exact,
  native_exact=native_exact,
  projected_strict=function(case) projected_original(case, 1e-7),
  projected_final=function(case) projected_original(case, 1e-6)
)

raw_results <- list()
comparison <- list()
for(case_name in names(cases)) {
  case <- cases[[case_name]]
  reference <- NULL
  for(method_name in names(methods)) {
    for(run in seq_len(repetitions)) {
      gc()
      result <- methods[[method_name]](case)
      raw_results[[length(raw_results) + 1L]] <- data.frame(
        case=case_name, D=length(case$documents), V=length(case$vocab),
        K=case$K, method=method_name, run=run,
        moment=result$moment, anchor=result$anchor,
        recovery=result$recovery, total=result$total)
      if(run == 1L && method_name == "legacy") reference <- result
      if(run == 1L) {
        comparison[[length(comparison) + 1L]] <- data.frame(
          case=case_name, method=method_name,
          anchor_overlap=length(intersect(reference$anchors, result$anchors)),
          max_beta_difference=if(method_name %in% c("legacy", "native_exact"))
            max(abs(reference$beta - result$beta)) else NA_real_)
      }
    }
  }
}

raw_results <- do.call(rbind, raw_results)
comparison <- do.call(rbind, comparison)
summary <- aggregate(
  cbind(moment, anchor, recovery, total) ~ case + D + V + K + method,
  raw_results, median)
summary <- merge(summary, comparison, by=c("case", "method"), all.x=TRUE)
reference_time <- setNames(
  summary$total[summary$method == "legacy"],
  summary$case[summary$method == "legacy"])
summary$speedup <- reference_time[summary$case]/summary$total
summary <- summary[order(summary$case, match(summary$method, names(methods))),]

print(summary, row.names=FALSE)
if(!is.null(output)) {
  saveRDS(list(summary=summary, raw=raw_results, comparison=comparison), output)
  cat("Saved benchmark results to", output, "\n")
}

