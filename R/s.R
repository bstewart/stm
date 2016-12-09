##
#Basis function generator
##
s <- function(x, df, ...) {
  if(class(x)=="Date") {
    warning("A Date object coerced to numeric. 
            Converting variable in advance will stop this warning in the future.
            Postprocessing tools may not work with dates.")
    x <- as.numeric(x)
  }
  nval <- length(unique(x))
  if(nval==1) stop("Smooth requested on covariate with only one value.")
  if(nval==2) {
    warning("Smooth requested on covariate with only two values. Converting to dichotomous coding.")
    return(as.numeric(as.factor(x)))
  }
  if(missing(df)) {
    df <- min(10, (nval-1))
  }
  return(bs(x, df,...))     
}


