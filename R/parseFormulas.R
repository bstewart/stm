#Function that parses formulas for plot.estimateEffect
#prep: output of estimateEffect
#cdata: matrix of control variables in the original form that will be used to
#be evaluated and then plotted.

parseFormulas <- function(prep, cdata){
  formula <- prep$formula
  #newdat <- rbind(prep$data, cdata)
  #elements of formula
  elements <- colnames(attr(terms(prep$formula), "factors"))
  #Identify splines, bsplines, etc
  special <- attr(terms(prep$formula, specials=c("s", "ns", "bs")),
  "special")
  #Replace old formula with new formula
  newform <- formula
  if(!is.null(special$s)){
    if(length(special$s)>1) stop("Only one spline supported")
    newform <- str_replace(as.character(newform)[2], "s\\(", "s2\\(")
    s2 <- function(x,...) predict(prep$modelframe[[elements[special$s]]],x)
  }
  if(!is.null(special$bs)){
    if(length(special$bs)>1) stop("Only one spline supported")
    newform <- str_replace(as.character(newform)[2], "bs\\(",
                           "bs2\\(")
    bs2 <- function(x,...) predict(prep$modelframe[[elements[special$bs]]],x)
  }
  if(!is.null(special$ns)){
    if(length(special$ns)>1) stop("Only one spline supported")
    newform <- str_replace(as.character(newform)[2], "ns\\(",
                           "ns2\\(")
    ns2 <- function(x,...) predict(prep$modelframe[[elements[special$ns]]],x)
  }
  if(is.null(special$ns) & is.null(special$s) & is.null(special$bs)){
    newform <- as.character(newform[2])
  }
  #Put it back in formula mode
  newform <- as.formula(paste("~", newform))
  #Make control matrix
  #newX <- model.matrix(newform, newdat)
  newX <- model.matrix(newform, cdata)
 # cdatout <- newX[(nrow(prep$data)+1):nrow(newX),]
  cdatout <- newX
  return(cdatout)
}
