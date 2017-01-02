#' Make a Design Matrix
#' 
#' Create a sparse model matrix which respects the basis functions
#' of the original data on which it was created. Primarily for internal
#' use but may be of some independent interest.
#' 
#' This functions is designed to be used in settings where we need
#' to make a prediction using a model matrix.  The practical challenge
#' here is ensuring that the representation of the data lines up
#' with the original representation.  This becomes challenging for
#' functions that produce a different representation depending on their
#' inputs. A simple conceptual example is factor variables.  If we run 
#' our original model using a factor with levels \code{c("A","B", "C")}
#' then when we try to make predictions for data having only levels
#' \code{c("A","C")} we need to adjust for the missing level.  Base
#' R functions like \code{\link{predict.lm}} in \pkg{stats} handle 
#' this gracefully and this function is essentially a version of 
#' \code{\link{predict.lm}} that only constructs the model matrix.
#' 
#' Beyond factors the key use case for this are basis functions like
#' splines.  For a function like this to work it must either depend only
#' on the observation it is transforming (e.g. \code{\link{log}}) or it must
#' have a generic for \code{\link{predict}} and \code{\link{makepredictcall}}.
#' The spline wrapper \code{\link{s}} has both and so should work.
#' 
#' When a function lacks these methods it will still produce a design matrix
#' but the values will be wrong.  To catch these settings we implement a quick
#' test when \code{test=TRUE} as it is by default.  To test we simply split
#' the original data in half and ensure that looking at each half separately produces
#' the same values as the complete original data.
#' 
#' @param formula the formula describing the design matrix.  Any responses will be deleted
#' @param origData the original dataset as a dataframe
#' @param newData a dataframe containing any of the variables in the formula.  This will provide
#' the data in the returned model matrix.
#' @param test when set to TRUE runs a test that the matrix was constructed correctly
#' see details for more.
#' @param sparse by default returns a sparse matrix using \code{\link[Matrix]{sparse.model.matrix}} in 
#' \pkg{Matrix}
#' @seealso \code{\link{fitNewDocuments}}
#' @examples 
#' foo <- data.frame(response=rnorm(30),
#'                   predictor=as.factor(rep(c("A","B","C"),10)), 
#'                   predictor2=rnorm(30))
#' foo.new <- data.frame(predictor=as.factor(c("A","C","C")), 
#'                       predictor2=foo$predictor2[1:3])
#' makeDesignMatrix(~predictor + s(predictor2), foo, foo.new)
#' @keywords internal
#' @export
makeDesignMatrix <- function(formula, 
                             origData, newData, 
                             test=TRUE, sparse=TRUE) {

  #this is inspired from the way predict.lm is structured 
  termobj <- terms(formula, data = origData)
  termobj <- stats::delete.response(termobj)
  mf <- model.frame(termobj, data = origData)
  mt <- attr(mf, "terms")
  mm <- model.matrix(mt, mf)
  contrasts <- attr(mm,"contrasts")
  xlevels <- stats::.getXlevels(mt,mf)
  newmf <- model.frame(mt, newData, xlev=xlevels)
  if(sparse) {
    X <- try(Matrix::sparse.model.matrix(mt, newmf, contrasts.arg = contrasts))
    if(class(X)=="try-error") X <- try(stats::model.matrix(mt,newmf, contrasts.arg=contrasts))
    if(class(X)=="try-error") stop("Error creating model matrix.")
    X <- Matrix::Matrix(X)
  } else {
    X <- try(model.matrix(mt, newmf, contrasts.arg = contrasts))
  }
  if(class(X)=="try-error") stop("Error creating model matrix.")
  
  if(test) {
    #test each half of the data
    halves <- list(seq(from=1, to=nrow(origData), by=2),
                   seq(from=2, to=nrow(origData), by=2))
    
    for(i in 1:2) {
      #create a test model frame and a test X
      tmf <- model.frame(mt, origData[halves[[i]],,drop=FALSE], xlev=xlevels)
      tX <- model.matrix(mt,tmf, contrasts.arg=contrasts)
      if(!isTRUE(all.equal(tX, mm[halves[[i]],,drop=FALSE], check.attributes=FALSE))) {
        stop("Model matrix constructed but testing failed.  This suggests you are using a function without
             a makepredictcall() and predict() generic.")
      }
      }
    }
  return(X)
  }

