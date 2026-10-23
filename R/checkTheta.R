check.theta.computed <- function(model, function.name) {
  if(!is.null(model$theta)) return(invisible(TRUE))
  stop(
    sprintf("%s requires model$theta, but theta was not computed.\n", function.name),
    "This happens for incremental inference on large corpora (or when ",
    "incremental.final.theta=FALSE). To compute document-topic proportions ",
    "at the fitted global parameters, run:\n",
    "  model <- updateTheta(model, documents)",
    call.=FALSE
  )
}
