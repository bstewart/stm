#plot.estimateEffect: Function to plot output from STM with respect to
#covariates
#x: output of estimateEffect
#covariate: covariate of interest
#model: model output
#topics: topics of interest
#method: plotting method, either pointestimate, difference, or continuous
# cov.value1: first value of covariate of interest if using the method
#"difference"
#cov.value2: second value of covariate of interest if using the method
#"difference"
#moderator: if using an interaction, variable that is the moderator of
#the covariate of interest
#moderator.value: if using an interaction, the value at which the moderator
#should be set
#npoints: for method continuous, the number of points that should be
#used for simulation on the continuous line
#nsims: the number of simulations used to create effect estimates
#the betas
#ci.level: confidence level for confidence intervals
#printlegend: whether or not to print a legend for method continuous
#labeltype: type of label
#n: number of words to print if labeltype is 'prob', 'frex', 'lift',
#or 'score
#frexw: frex weight if labeltype is 'frex'
#add: make a new plot (F) or add to existing plot (T)
#linecol: vector of line colors (should be the length of the number of
#topics, if method="continous"
#verbose.labels: adds covariate information to the labels
#width: width of labels on the plot
#family: font family
#custom.labels: vector of custom labels if labeltype="custom"

#' Plot effect of covariates on topics
#' 
#' Plots the effect of a covariate on a set of topics selected by the user.
#' Different effect types available depending on type of covariate. Before
#' running this, the user should run a function to simulate the necessary
#' confidence intervals.  See \code{\link{estimateEffect}}.
#' 
#' 
#' @param x Output of estimateEffect, which calculates simulated betas for
#' plotting.
#' @param covariate String of the name of the main covariate of interest. Must
#' be enclosed in quotes.  All other covariates within the formula specified in
#' estimateEffect will be kept at their median.
#' @param model Model output, only necessary if labeltype is "prob", "frex",
#' "score", or "lift".  Models with more than one spline cannot be used for
#' plot.estimateEffect.
#' @param topics Topics to plot.
#' @param method Method used for plotting.  "pointestimate" estimates mean
#' topic proportions for each value of the covariate.  "difference" estimates
#' the mean difference in topic proportions for two different values of the
#' covariate (cov.value1 and cov.value2 must be specified).  "continuous"
#' estimates how topic proportions vary over the support of a continuous
#' covariate.
#' @param cov.value1 For method "difference", the value or set of values of
#' interest at which to set the covariate. In the case of calculating a
#' treatment/control contrast, set the treatment to cov.value1.
#' @param cov.value2 For method "difference", the value or set of values which
#' will be set as the comparison group.  cov.value1 and cov.value2 must be
#' vectors of the same length.
#' @param moderator When two terms are interacted and one variable in the
#' interaction is the covariate of interest, the user can specify the value of
#' the interaction with moderator.value, and the name of the moderator with
#' moderator.
#' @param moderator.value When two terms are interacted and one variable in the
#' interaction is the covariate of interest, the user can specify the value of
#' the interaction term.
#' @param npoints Number of unique points to use for simulation along the
#' support of a continuous covariate.  For method "continuous" only.
#' @param nsims Number of simulations for estimation.
#' @param n Number of words to print if "prob", "score", "lift", or "frex" is
#' chosen.
#' @param ci.level Confidence level for confidence intervals.
#' @param frexw If "frex" labeltype is used, this will be the frex weight.
#' @param add Logical parameter for whether the line should be added to the
#' plot, or a new plot should be drawn.
#' @param linecol For continuous covariates only.  A vector that specifies the
#' colors of the lines within the plot.  If NULL, then colors will be randomly
#' generated.
#' @param verbose.labels For method "difference" -- verboselabels will specify
#' the comparison covariate values of the covariate on the plot.
#' @param xlim Vector of x axis minimum and maximum values.
#' @param ylim Vector of y axis minimum and maximum values.
#' @param main Character string that is plot title.
#' @param printlegend Whether to plot a topic legend in the case of a
#' continuous covariate.
#' @param labeltype Determines the labeltype for the topics.  The default is
#' "number" which prints the topic number.  Other options are "prob", which
#' prints the highest probability words, "score", "lift", and "frex", from
#' labeltopics (see labeltopics() for more details).  The user can also select
#' "custom" for custom labels, which should be inputted under custom.labels.
#' Labels appear in the legend for continous covariates.
#' @param xlab Character string that is x axis title.
#' @param ylab Character string that is y axis title.
#' @param width Number that specifies width of the character string.  Smaller
#' numbers will have smaller-width labels.  Default is 25.
#' @param custom.labels A vector of custom labels if labeltype is equal to
#' "custom".
#' @param family Font family.
#' @param ...  Other plotting parameters
#' @return Values returned invisibly will depend on the method
#'  
#' For pointestimate: 
#' \item{uvals}{Values of the covariate at which means and ci's
#' were evaluated.} 
#' \item{topics}{Topics for which means and ci's were
#' evaluated.} 
#' \item{means}{For each topic, means for each unique value.}
#' \item{cis}{For each topic, confidence intervals for each unique value.}
#' \item{labels}{Labels for each topic and unique value.}
#' 
#' For difference: 
#' \item{topics}{Topics for which difference in means and ci's
#' were evaluated} 
#' \item{means}{For each topic, difference in means.}
#' \item{cis}{For each topic, confidence intervals for difference in means.}
#' \item{labels}{Labels for each topic.}
#' 
#' For continuous: 
#' \item{x}{Individual values of the covariate at which means
#' and ci's were evaluated.} 
#' \item{topics}{Topics for which means and ci's
#' were evaluated} 
#' \item{means }{For each topic and each x, means.} 
#' \item{cis}{For each topic and each x, confidence intervals for difference in means.}
#' \item{labels}{Labels for each topic.}
#' 
#' @examples
#' 
#' \dontrun{
#' 
#' prep <- estimateEffect(1:3 ~ treatment, gadarianFit, gadarian)
#' plot(prep, "treatment", model=gadarianFit,
#' method="pointestimate")
#' plot(prep, "treatment", model=gadarianFit,
#' method="difference",cov.value1=1,cov.value2=0)
#' 
#' #If the covariate were a binary factor, 
#' #the factor labels can be used to  
#' #specify the values of cov.value1 (e.g., cov.value1="treat"). 
#' 
#' # String variables must be turned to factors prior to plotting. 
#' #If you see this error, Error in rep.int(c(1, numeric(n)), n - 1L) : 
#' # invalid 'times' value, then you likely have not done this.
#' 
#' #Example of binary times binary interaction
#' gadarian$binaryvar <- sample(c(0,1), nrow(gadarian), replace=T)
#' temp <- textProcessor(gadarian$open.ended.response,metadata=gadarian)
#' out <- prepDocuments(temp$documents, temp$vocab, temp$meta)
#' stm1 <- stm(out$documents, out$vocab, 3, prevalence=~treatment*binaryvar,
#'  data=gadarian)
#' prep <- estimateEffect(c(2) ~ treatment*binaryvar, stmobj=stm1,
#' metadata=gadarian)
#' 
#' par(mfrow=c(1,2))
#' plot(prep, "treatment", method="pointestimate",
#' cov.value1=1, cov.value2=0, xlim=c(-1,1), moderator="binaryvar", moderator.value=1)
#' plot(prep, "treatment", method="pointestimate",
#' cov.value1=1, cov.value2=0, xlim=c(-1,1), moderator="binaryvar",
#' moderator.value=0)
#' }
#' @export 
plot.estimateEffect <- function(x, covariate, model=NULL,
                                topics=x$topics,
                                method=c("pointestimate", "difference","continuous"),
                                cov.value1=NULL, cov.value2=NULL,
                                moderator=NULL, moderator.value=NULL,
                                npoints=100, nsims=100, ci.level=.95,
                                xlim=NULL, ylim=NULL, xlab="",ylab=NULL,
                                main="", printlegend=T,
                                labeltype="numbers", n=7, frexw=.5,
                                add=F, linecol=NULL, width=25,
                                verbose.labels=T, family=NULL,
                                custom.labels=NULL,...){
  
  method <- match.arg(method)
  if(method=="difference" && (is.null(cov.value1) | is.null(cov.value2))) {
    stop("For method='difference' both cov.value1 and cov.value2 must be specified.")
  }
  #Produce cdata (data in original form) and
  #cmatrix (data in design matrix form)
  cthis <- produce_cmatrix(prep=x, covariate=covariate, method=method,
                           cov.value1=cov.value1,
                           cov.value2=cov.value2, npoints=npoints,
                           moderator=moderator, moderator.value=moderator.value)
  cdata <- cthis$cdata
  cmat <- cthis$cmatrix

  #Simulate betas
  simbetas <- simBetas(x$parameters, nsims=nsims)

  #Find offset for confidence level
  offset <- (1-ci.level)/2

  
  #Plot for each method
  if(method=="continuous"){
    toreturn <- plotContinuous(prep=x,covariate=covariate,topics=topics, cdata=cdata, cmat=cmat, simbetas=simbetas,
                   offset=offset,xlab=xlab, ylab=ylab, main=main,
                   xlim=xlim, ylim=ylim, linecol=linecol, add=add,
                   labeltype=labeltype,n=n,custom.labels=custom.labels,model=model,frexw=frexw,printlegend=printlegend,...)
    return(invisible(toreturn))
  }
  if(method=="pointestimate"){
    toreturn <- plotPointEstimate(prep=x,covariate=covariate,topics=topics, cdata=cdata, cmat=cmat, simbetas=simbetas,
                      offset=offset,xlab=xlab, ylab=ylab, main=main,
                      xlim=xlim, ylim=ylim, linecol=linecol, add=add,
                      labeltype=labeltype,n=n,
                                  custom.labels=custom.labels,model=model,frexw=frexw,width=width,
                                  verbose.labels=verbose.labels,...)
    return(invisible(toreturn))
  }
  if(method=="difference"){
    if(missing(cov.value1)) stop("Missing a value for cov.value1. See documentation.")
    if(missing(cov.value2)) stop("Missing a value for cov.value2. See documentation.")
    toreturn <- plotDifference(prep=x,covariate=covariate,topics=topics, cdata=cdata, cmat=cmat, simbetas=simbetas,
                   offset=offset,xlab=xlab, ylab=ylab, main=main,
                   xlim=xlim, ylim=ylim, linecol=linecol, add=add,
                   labeltype=labeltype,n=n,
                   custom.labels=custom.labels,model=model,frexw=frexw,width=width,
                   cov.value1=cov.value1,
                               cov.value2=cov.value2,verbose.labels=verbose.labels,...)
    return(invisible(toreturn))
  }
}