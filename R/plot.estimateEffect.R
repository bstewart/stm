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


plot.estimateEffect <- function(x, covariate, model=NULL,
                                topics=x$topics,
                                method="pointestimate",
                                cov.value1=NULL, cov.value2=NULL,
                                moderator=NULL, moderator.value=NULL,
                                npoints=100, nsims=100, ci.level=.95,
                                xlim=NULL, ylim=NULL, xlab="",ylab=NULL,
                                main="", printlegend=T,
                                labeltype="numbers", n=7, frexw=.5,
                                add=F, linecol=NULL, width=25,
                                verbose.labels=T, family=NULL,
                                custom.labels=NULL,...){

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
