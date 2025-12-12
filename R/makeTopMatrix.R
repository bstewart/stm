#' Convert prevalence covariates to design matrix
#'
#' Internal function to process prevalence formulas or matrices into
#' design matrices for STM models. Handles both sparse and dense matrices
#' and ensures proper intercept column.
#'
#' @param x Either a formula or matrix of prevalence covariates
#' @param data Data frame containing variables referenced in formula
#' @return Design matrix (sparse or dense depending on sparsity)
#' @keywords internal
#' @noRd
makeTopMatrix <- function(x, data=NULL) {
  #is it a formula?
  if(inherits(x,"formula")) {
    termobj <- terms(x, data=data)
    if(attr(termobj, "response")==1) stop("Response variables should not be included in prevalence formula.")
    xmat <- try(Matrix::sparse.model.matrix(termobj,data=data),silent=TRUE)
    if(inherits(xmat,"try-error")) {
      xmat <- try(stats::model.matrix(termobj, data=data), silent=TRUE)
      if(inherits(xmat,"try-error")) {
               stop("Error creating model matrix.
               This could be caused by many things including
               explicit calls to a namespace within the formula.
               Try a simpler formula.")
      }
      xmat <- Matrix::Matrix(xmat)
    }
    propSparse <- 1 - Matrix::nnzero(xmat)/length(xmat)
    #if its less than 50% sparse or there are fewer than 50 columns, just convert to a standard matrix
    if(propSparse < .5 | ncol(xmat) < 50) {
      xmat <- as.matrix(xmat)
    }
    return(xmat)
  }
  if(is.matrix(x)) {
    #Does it have an intercept in first column?
    if(isTRUE(all.equal(x[,1],rep(1,nrow(x))))) return(Matrix::Matrix(x))
    else return(cbind(1,Matrix::Matrix(x)))
  }
}
