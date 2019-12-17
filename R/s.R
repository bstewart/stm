#' Make a B-spline Basis Function
#' 
#' This is a simple wrapper around the \code{\link[splines]{bs}} function in
#' the splines package.  It will default to a spline with 10 degrees of
#' freedom.
#' 
#' This is a simple wrapper written as users may find it easier to simply type
#' \code{s} rather than selecting parameters for a spline. We also include
#' \code{predict} and \code{makepredictcall} generic functions for the class
#' so it will work in settings where \code{\link{predict}} is called.
#' 
#' @param x The predictor value.
#' @param df Degrees of freedom.  Defaults to the minimum of 10 or one minus
#' the number of unique values in x.
#' @param \dots Arguments passed to the \code{\link[splines]{bs}} function.
#' @return A predictor matrix of the basis functions.
#' @seealso \code{\link[splines]{bs}} \code{\link[splines]{ns}}
#' @export
s <- function(x, df, ...) {
  if(inherits(x,"Date")) {
    warning("A Date object coerced to numeric. 
            Converting variable in advance will stop this warning in the future.
            Postprocessing tools may not work with dates.")
    x <- as.numeric(x)
  }
  
  nval <- length(unique(x))
  if(missing(df)) {
    df <- min(10, (nval-1))
  }
  obj <- splines::bs(x, df,...)
  attr(obj, "class") <- c("s", attr(obj, "class")) #we need this to ensure that our predict generics trigger
  return(obj)
}

#' @export
#' @keywords internal
predict.s <- function (object, newx, ...) 
{
  if (missing(newx)) 
    return(object)
  a <- c(list(x = newx), attributes(object)[c("degree", "knots", 
                                              "Boundary.knots", "intercept")])
  do.call("splines::bs", a)
}

#' @export
#' @importFrom stats makepredictcall
#' @keywords internal
makepredictcall.s <- function (var, call) 
{
  #if (as.character(call)[1L] != "bs") 
  #  return(call)
  at <- attributes(var)[c("degree", "knots", "Boundary.knots", 
                          "intercept")]
  xxx <- call[1L:2]
  xxx[names(at)] <- at
  xxx
}
