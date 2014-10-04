
permutationTest <- function(formula, stmobj, treatment, 
                         nruns=100, 
                         documents, vocab, data, seed=NULL,
                         stmverbose=TRUE,uncertainty="Global") {
  if(!requireNamespace("clue", quietly=TRUE)) stop("Install the clue package to use this function")
  
  settings <- stmobj$settings
  
  if(!is.data.frame(data)) stop("data object must be a data.frame containing treatment variable.")
  if(!(treatment %in% colnames(data))) stop("treatment variable must be in the data.frame")
  if(!all(data[,treatment]%in%c(0,1))) stop("treatment variable must be binary and coded 0 and 1.")
  prob <- mean(data[,treatment])

  #use the same seed behavior as stm()
  if(is.null(seed)) {
    seed <- floor(runif(1)*1e7) 
    set.seed(seed)
  } else {
    set.seed(seed)
  }
  
  #manipulate the settings object a bit
  settings$verbose <- stmverbose
  settings$keepHistory <- FALSE
  
  #Save the model beta for alignment
  betaref <- exp(stmobj$beta$logbeta[[1]])
  
  #a small internal function to automate the calculation of the quantile of effects
  qeffects <- function(formula, mod, data, uncertainty) {
    #estimate the model effect
    prep <- estimateEffect(formula, stmobj=mod, metadata=data, uncertainty=uncertainty)
    betas <- lapply(prep$parameters, function (x) lapply(x, function (y) rmvnorm(100, y$est, y$vcov)))
    simbetas <- lapply(betas, do.call, what=rbind)
    #which column is it in the parameter matrix?
    parcol <- which(colnames(settings$covariates$X)==treatment)
    par <- lapply(simbetas, function(x) quantile(x[,parcol], c(.025, .5, .975)))
    par
  }
  
  #start with the reference model effect
  cat("Calculating effects for reference model \n")
  ref.effect <- qeffects(formula, stmobj, data, uncertainty)
  
  #now do the additional runs
  tosave <- list()
  for(i in 1:(nruns-1)) {
    #set the seed and randomize the treatment
    settings$seed <- floor(runif(1)*1e7)
    data[,treatment] <- rbinom(n=nrow(data), size=1, prob=prob)
    #now copy the new treatment value into the data
    termobj <- terms(formula, data=data)
    if(attr(termobj, "response")==1) stop("Response variables should not be included in prevalence formula.")
    settings$covariates$X <- model.matrix(termobj,data=data)
    #run the model
    cat(sprintf("Running model %i of %i \n", (i+1), (nruns)))
    mod <- stm.control(documents, vocab, settings=settings)
    par <- qeffects(formula, mod, data, uncertainty)
    
    betamod <- exp(mod$beta$logbeta[[1]])
    align <- clue::solve_LSAP(betaref%*%t(betamod), maximum=TRUE) 
    tosave[[i]] <- par[align]
  }
  out <- list(ref=ref.effect, permute=tosave, variable=treatment, seed=seed)
  class(out) <- "STMpermute"
  return(out)
}

plot.STMpermute <- function(x, topic,
                            type=c("match", "largest"),
                            xlim=NULL, ylim=NULL, ...) {
  permuteobj <- x
  type <- match.arg(type)
  if(length(topic)>1) stop("topic must be a single integer.")
  ref <- as.numeric(permuteobj$ref[[topic]])
  placebo <- permuteobj$permute
  #unpack all the placebo elements
  if(type=="match") {
    placebo <- do.call(rbind,unlist(lapply(placebo, function(x) x[topic]), recursive=FALSE))
  } else {
    maxfunc <- function(x) {
      x[[which.max(sign(ref[2])*do.call(rbind,x)[,2])]]
    }
    placebo <- do.call(rbind,lapply(placebo, maxfunc))
  }
  toplot <- as.matrix(rbind(placebo, ref))
  if(is.null(xlim)) xlim <- c(min(toplot), max(toplot))
  if(is.null(ylim)) ylim <- c(0, nrow(toplot)+1)

  plot(toplot[,2], 1:nrow(toplot), xlim=xlim, ylim=ylim, ...)
  n1 <- nrow(toplot)-1
  segments(toplot[1:n1,1],1:n1, toplot[1:n1,3], 1:n1)  
  n <- n1+1
  segments(toplot[n,1], n, toplot[n,3], n, col="red")
  abline(v=0, lty=2)
}