#' Permutation test of a binary covariate.
#' 
#' Run a permutation test where a binary treatment variable is randomly
#' permuted and topic model is reestimated.
#' 
#' This function takes a single binary covariate and runs a permutation test
#' where, rather than using the true assignment, the covariate is randomly
#' drawn with probability equal to its empirical probability in the data. After
#' each shuffle of the covariate the same STM model is estimated at different
#' starting values using the same initialization procedure as the original
#' model, and the effect of the covariate across topics is calculated.
#' 
#' Next the function records two quantities of interest across this set of
#' "runs" of the model. The first records the absolute maximum effect of the
#' permuted covariate across all topics.
#' 
#' The second records the effect of the (permuted) covariate on the topic in
#' each additional stm run which is estimated to be the topic closest to the
#' topic of interest (specified in \code{\link{plot.STMpermute}}) from the
#' original stm model. Uncertainty can be calculated using the standard options
#' in \code{\link{estimateEffect}}.
#' 
#' @param formula A formula for the prevalence component of the \code{stm}
#' model and the \code{estimateEffect} call.  This formula must contain at
#' least one binary covariate (specified through the argument \code{treatment})
#' but it can contain other terms as well.  If the binary covariate is
#' interacted with additional variables the estimated quantity of interest is
#' the effect when those additional variables are set to 0.
#' @param stmobj Model output from a single run of \code{stm} which contains
#' the reference effect.
#' @param treatment A character string containing treatment id as used in the
#' formula of the stmobj.  This is the variable which is randomly permuted.
#' @param nruns Number of total models to fit (including the original model).
#' @param documents The documents used in the stmobj model.
#' @param vocab The vocab used in the stmobj model.
#' @param data The data used in the stmobj model.
#' @param seed Optionally a seed with which to replicate the result.  As in
#' \code{\link{stm}} the seed is automatically saved and returned as part of
#' the object.  Passing the seed here will replicate the previous run.
#' @param stmverbose Should the stm model be run with \code{verbose=TRUE}.
#' Turning this to \code{FALSE} will surpress only the model specific printing.
#' An update on which model is being run will still print to the screen.
#' @param uncertainty Which procedure should be used to approximate the
#' measurement uncertainty in the topic proportions.  See details for more
#' information.  Defaults to the Global approximation.
#' @return \item{ref}{A list of K elements containing the quantiles of the
#' estimated effect for the reference model.} \item{permute}{A list where each
#' element is an aligned model parameter summary} \item{variable}{The variable
#' id that was permuted.} \item{seed}{The seed for the stm model.}
#' @seealso \code{\link{plot.STMpermute}}
#' @examples
#' 
#' \dontrun{
#' temp<-textProcessor(documents=gadarian$open.ended.response,metadata=gadarian)
#' out <- prepDocuments(temp$documents, temp$vocab, temp$meta)
#' documents <- out$documents
#' vocab <- out$vocab
#' meta <- out$meta
#' set.seed(02138)
#' mod.out <- stm(documents, vocab, 3, prevalence=~treatment + s(pid_rep), data=meta)
#' summary(mod.out)
#' prep <- estimateEffect(1:3 ~ treatment + s(pid_rep), mod.out, meta)
#' plot(prep, "treatment", model=mod.out,
#'      method="difference",cov.value1=1,cov.value2=0)
#' test <- permutationTest(formula=~ treatment + s(pid_rep), stmobj=mod.out, 
#'                         treatment="treatment", nruns=25, documents=documents,
#'                         vocab=vocab,data=meta, stmverbose=FALSE)
#' plot(test,2, xlab="Effect", ylab="Model Index", main="Topic 2 Placebo Test")
#' }
#' @export
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
    mod <- stm.control(documents, vocab, settings=settings, model=NULL)
    par <- qeffects(formula, mod, data, uncertainty)
    
    betamod <- exp(mod$beta$logbeta[[1]])
    align <- clue::solve_LSAP(betaref%*%t(betamod), maximum=TRUE) 
    tosave[[i]] <- par[align]
  }
  out <- list(ref=ref.effect, permute=tosave, variable=treatment, seed=seed)
  class(out) <- "STMpermute"
  return(out)
}

#' Plot an STM permutation test.
#' 
#' Plots the results of a permutation test run using
#' \code{\link{permutationTest}}.
#' 
#' This function plots the output of \code{\link{permutationTest}} by stacking
#' horizontal confidence intervals for the effects of the permuted variable.
#' In choosing the topic in the permuted runs of stm to plot the effect for,
#' two methods are available, "match" and "largest". The former uses Kuhn's
#' (1955) Hungarian method to align the topics, and then uses the model's best
#' match of the reference topic.  The latter uses the topic which has the
#' expected effect size in the direction of the reference model effect; thus,
#' we would expect this method to be quite conservative.
#' 
#' @param x Object from the output of \code{\link{permutationTest}}.
#' @param topic Integer indicating which topic to plot.
#' @param type Character string indicating what topic comparison to use.
#' "match" uses the Hungarian aligned method and "largest" uses the largest
#' mean in direction of reference topic.
#' @param xlim Range of the X-axis.
#' @param ylim Range of the Y-axis.
#' @param ...  Other parameters which may be passed to plot.
#' @seealso \code{\link{permutationTest}}
#' @examples
#' 
#' \dontrun{
#' 
#' temp<-textProcessor(documents=gadarian$open.ended.response,metadata=gadarian)
#' out <- prepDocuments(temp$documents, temp$vocab, temp$meta)
#' documents <- out$documents
#' vocab <- out$vocab
#' meta <- out$meta
#' set.seed(02138)
#' mod.out <- stm(documents, vocab, 3, prevalence=~treatment + s(pid_rep), data=meta)
#' summary(mod.out)
#' prep <- estimateEffect(1:3 ~ treatment + s(pid_rep), mod.out, meta)
#' plot(prep, "treatment", model=mod.out,
#'      method="difference",cov.value1=1,cov.value2=0)
#' test <- permutationTest(formula=~ treatment + s(pid_rep), stmobj=mod.out, 
#'                         treatment="treatment", nruns=25, documents=documents,
#'                         vocab=vocab,data=meta, stmverbose=FALSE)
#' plot(test,2, xlab="Effect", ylab="Model Index", main="Topic 2 Placebo Test")
#' }
#' @export
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