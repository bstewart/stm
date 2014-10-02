plot.estimateEffect <- function(x, covariate, model=NULL, topics=x$topics,
                          method="pointestimate",
                          cov.value1=NULL, cov.value2=NULL,int.value=NULL, 
                          npoints=100, nsims=100, ci.level=.95, 
                          xlim=NULL, ylim=NULL, ylab=NULL, 
                          main="",printlegend=TRUE, 
                          labeltype="numbers",n=7,frexw=.5, 
                          xlab="", add=FALSE, linecol=NULL, 
                          width=25, verbose.labels=TRUE, 
                          family=NULL, custom.labels=NULL,...){


##################
##Preliminaries###
##################
object <- x

#Grab the model matrix, number of params in the model
xmat <- model.matrix(object$formula, data=object$data)
nparams <- ncol(xmat)

#Are there any interactions?
terms <- attr(object$modelframe, "terms")
factors <- attr(terms, "factors")
fnames <- colnames(factors)
parse <- grepl(":", fnames, fixed = TRUE)

#Are there any factors?
facts <- attr(xmat, "contrasts")
factornames <- names(facts)

#Error checking -- did they put quotes around covariate
if(!is.character(covariate)){
  stop("Please enter the name of the covariate of interest as a string.")
}

#Interaction error checking, find out whether the covariate of interest is in the interaction.
if(any(parse)) {
    intterms <- which(parse)
    if(sum(parse)>1){
      stop("Plot only supports one interaction term.")
    }
    intlabs <- c(str_split(fnames[intterms], ":")[[1]][1], str_split(fnames[intterms], ":")[[1]][2])
    if(length(unique(object$data[,intlabs[1]]))>2 & length(unique(object$data[,intlabs[2]]))>2){
      stop("Plot only supports interaction terms where one of the covariates has two values.")
    }
    #if(!is.numeric(object$data[,intlabs[1]]) & !is.numeric(object$data[,intlabs[2]])){
    #  stop("Plot only supports interaction between two numeric variables.")
    #}
    if(covariate%in%intlabs){
      covint <- which(intlabs==covariate)
      cint <- which(intlabs!=covariate)
    }else{
      covint <- 0
    }
  }

#Find the splines
classes <- attr(attr(object$modelframe, "terms"), "dataClasses")
splines <- sapply(classes, function(x) grep("nmatrix", x, fixed=T, value=T))
splinesvars <- names(unlist(splines))

######################
##Simulate the betas##
#####################
betas <- lapply(object$parameters, function (x) lapply(x, function (y) rmvnorm(nsims, y$est, y$vcov)))
simbetas <- list()
for(i in 1:length(betas)){
  temp <- matrix(ncol=nparams)
  for(j in 1:length(betas[[i]])){
    temp <- rbind(temp, betas[[i]][[j]])
  }
  simbetas[[i]] <- temp[-1,]
}

#################################
##Calculate the control vector##
#################################
labels1 <- object$varlist
term.labels <- attr(terms(object$modelframe), "term.labels")
#NOTE: this will only be right for the median, if we change to the mean, we'll have to change this.
cvector <- apply(xmat,2,median)
matnames <- attr(xmat, "assign")

#Identify which entry in the control vector is the covariate of interest
covid <- which(matnames==which(labels1==covariate))

#Get control vector if there are splines in the function
if(!is.null(splinesvars)){
  for(i in 1:length(splinesvars)){
    splineid <- which(matnames==which(term.labels==splinesvars[i]))
    cvector[splineid] <- predict(object$modelframe[[splinesvars[i]]], median(object$data[,c(labels1[which(term.labels==splinesvars[i])])]))
    if(sum(covid==splineid)==length(covid) & method=="continuous"){
      method <- "spline"
      covspline <- splinesvars[i]
    }
  }
}

#Get the control vector if there are interactions in the function
if(any(parse)){
  cvector[intterms+1] <- median(xmat[,which(matnames==which(labels1==intlabs[1]))])*median(xmat[,which(matnames==which(labels1==intlabs[2]))])
}

##########
#Labeling#
##########

#Grab the labels
if(labeltype=="numbers"){
  labels=paste("Topic", object$topics)
}
if(labeltype=="prob"){
  if(is.null(model)){
    stop("Model output is needed in order to use words as labels.  Please enter the model output for the argument 'model'.")
  }
  else{
    labels = paste(apply(labelTopics(model, topics=object$topics, n=n)$prob[object$topics,],1,paste, collapse=","))
  }
}
if(labeltype=="lift"){
  if(is.null(model)){
    stop("Model output is needed in order to use words as labels.  Please enter the model output for the argument 'model'.")
  }
  else{
    labels = paste(apply(labelTopics(model, topics=object$topics, n=n)$lift[object$topics,],1,paste, collapse=","))
  }
}
if(labeltype=="frex"){
  if(is.null(model)){
    stop("Model output is needed in order to use words as labels.  Please enter the model output for the argument 'model'.")
  }
  else{
    labels = paste(apply(labelTopics(model, topics=object$topics, n=n,frexweight=frexw)$frex[object$topics,],1,paste, collapse=","))
  }
}
if(labeltype=="score"){
  if(is.null(model)){
    stop("Model output is needed in order to use words as labels.  Please enter the model output for the argument 'model'.")
  }
  else{
    labels = paste(apply(labelTopics(model, topics=object$topics, n=n)$score[object$topics,],1,paste, collapse=","))
  }
}
if(labeltype=="custom"){
  labels = rev(custom.labels)
}

###################################
#Determine confidence level offset#
###################################
offset <- (1-ci.level)/2

########################
#Point estimate method#
######################
if(method=="pointestimate"){
  #If the covariate is a factor
  if(length(covid)>1 & covariate%in%factornames){
    #Estimate simulated values to plot
    uvals <- xmat[,covid][!duplicated(xmat[,covid]),]
    toplot <- list()
    for(j in 1:nrow(uvals)){
      toplot[[j]] <- list()
      for(i in 1:length(object$topics)){
        cvector[covid] <- uvals[j,]
        if(any(parse)){
              #Interaction terms in xmat
          intterms2 <- str_detect(colnames(xmat), ":")
                                        #Interaction levels
          interactlevels <- unique(xmat[,which(matnames==which(labels1==intlabs[cint]))])  

          if(covint>0 & is.null(int.value)){
            cvector[which(matnames==which(labels1==intlabs[cint]))] <- median(xmat[,which(matnames==which(labels1==intlabs[cint]))])
            cvector[intterms2] <- uvals[j,]*median(xmat[,which(matnames==which(labels1==intlabs[cint]))])
          }
          if(covint>0 & !is.null(int.value)){
            cvector[which(matnames==which(labels1==intlabs[cint]))] <- int.value
            cvector[intterms2] <- uvals[j,]*int.value
          }
        }
        print(cvector)
        sims <- simbetas[[i]]%*%cvector
        toplot[[j]][[i]] <- list(mean=mean(sims), cis=quantile(sims, c(offset, 1-offset)))
      }
    }

    #Add string to labels to specify which covariate level is being plotted.
    newlabels <- list()
    
    factorlevs <- contr.treatment(levels(as.factor(object$data[,covariate])))
    for(i in 1:nrow(uvals)){
        newlabels[[i]] <- vector(length=length(object$topics))
        for(k in 1:length(object$topics)){
          if(verbose.labels){
            lab <- which(apply(factorlevs,1,function(x) sum(x==uvals[i,])==length(x)))
            newlabels[[i]][k] <- paste(labels[k], " (Covariate Level: ", rownames(factorlevs)[lab], ")",sep="")
          } else {
            newlabels[[i]][k] <- paste(labels[k])
          }
        }
      }
    
    #Plot estimates
    topicid <- rev(which(object$topics%in%topics))
    if (is.null(xlim) & length(labeltype)==1) if(labeltype!="numbers") xlim <- c(min(unlist(toplot)) - 1.5*abs(min(unlist(toplot))), max(unlist(toplot)))
    if (is.null(xlim) & length(labeltype)==1) if(labeltype=="numbers") xlim <- c(min(unlist(toplot)) - .2*abs(min(unlist(toplot))), max(unlist(toplot)))
    if (is.null(xlim) & length(labeltype)!=1) xlim <- c(min(unlist(toplot)) - .5*abs(min(unlist(toplot))), max(unlist(toplot)))
    if (is.null(ylim)) ylim <- c(0,length(topicid)*nrow(uvals)+1)
    if (is.null(ylab)) ylab <- ""
    if (add==F) plot(0,0, col="white", xlim=xlim, ylim=ylim, yaxt="n", ylab=ylab, xlab=xlab, main=main,...)
    lines(c(0,0), c(0, length(topicid)*nrow(uvals)+2), lty=2)
    it = 1
    for(j in 1:length(toplot)){
        for(i in 1:length(topicid)){
          points(toplot[[j]][[topicid[i]]]$mean, (j-1)*length(topicid)+i, pch=16)
          lines(c(toplot[[j]][[topicid[i]]]$cis[1], toplot[[j]][[topicid[i]]]$cis[2]), c((j-1)*length(topicid)+i,(j-1)*length(topicid)+i))
          if(labeltype!="custom") axis(2, at=(j-1)*length(topicid)+i, labels=str_wrap(newlabels[[j]][topicid[i]],width=width), las=1, cex=.25, tick=F, pos=toplot[[j]][[topicid[i]]]$cis[1])
          if(labeltype=="custom"){
            axis(2, at=(j-1)*length(topicid)+i, labels=str_wrap(labels[it],width=width), las=1, cex=.25, tick=F, pos=toplot[[j]][[topicid[i]]]$cis[1])
            it = it + 1
          }
        }
      }
  }else{
    #If the covariate is not a factor
    #Simulate quantities of interest
      uvals <- unique(xmat[,covid])
      if(length(uvals)==length(xmat[,covid])) stop("Too many levels to use point estimate method")
      toplot <- list()
      for(j in 1:length(uvals)){
        toplot[[j]] <- list()
        for(i in 1:length(object$topics)){
          cvector[covid] <- uvals[j]
          if(any(parse)){
            intcontrol <- which(colnames(xmat)==intlabs[cint])
            if(covint>0 & is.null(int.value)){
              cvector[intterms+1] <- uvals[j]*median(xmat[,which(matnames==which(labels1==intlabs[cint]))])
            }
            if(covint>0 & !is.null(int.value)){
              cvector[intcontrol] <- int.value
              cvector[intterms+1] <- uvals[j]*int.value
            }
          }
          sims <- simbetas[[i]]%*%cvector
          toplot[[j]][[i]] <- list(mean=mean(sims), cis=quantile(sims, c(offset, 1-offset)))
        }
      }
      #Add specific covariate levels to labels
      newlabels <- list()
      if(covariate%in%factornames){
        factorlevs <- contr.treatment(levels(as.factor(object$data[,covariate])))
      }
      for(i in 1:length(uvals)){
        newlabels[[i]] <- vector(length=length(object$topics))
        for(k in 1:length(object$topics)){
          if(verbose.labels){
            if(covariate%in%factornames){
              newlabels[[i]][k] <- paste(labels[k], " (Covariate Level: ", rownames(factorlevs)[factorlevs[,1]==uvals[i]], ")",sep="")
            }else{
              newlabels[[i]][k] <- paste(labels[k], " (Covariate Level: ", uvals[i], ")",sep="")
            }
          } else {
            newlabels[[i]][k] <- paste(labels[k])
          }
        }
      }

      #Plot estimates
      topicid <- rev(which(object$topics%in%topics))
      if (is.null(xlim)) if(labeltype!="numbers") xlim <- c(min(unlist(toplot)) - 1.5*abs(min(unlist(toplot))), max(unlist(toplot)))
      if (is.null(xlim)) if(labeltype=="numbers") xlim <- c(min(unlist(toplot)) - .2*abs(min(unlist(toplot))), max(unlist(toplot)))
      if (is.null(ylim)) ylim <- c(0,length(topicid)*length(uvals)+1)
      if (is.null(ylab)) ylab <- ""
      if (add==F) plot(0,0, col="white", xlim=xlim, ylim=ylim, yaxt="n", ylab=ylab, xlab=xlab, main=main,...)
      lines(c(0,0), c(0, length(topicid)*length(uvals)+2), lty=2)
      it = 1
      for(j in 1:length(toplot)){
        for(i in 1:length(topicid)){
          points(toplot[[j]][[topicid[i]]]$mean, (j-1)*length(topicid)+i, pch=16)
          lines(c(toplot[[j]][[topicid[i]]]$cis[1], toplot[[j]][[topicid[i]]]$cis[2]), c((j-1)*length(topicid)+i,(j-1)*length(topicid)+i))
          if(labeltype!="custom") axis(2, at=(j-1)*length(topicid)+i, labels=str_wrap(newlabels[[j]][topicid[i]],width=width), las=1, cex=.25, tick=F, pos=toplot[[j]][[topicid[i]]]$cis[1])
          if(labeltype=="custom"){
            axis(2, at=(j-1)*length(topicid)+i, labels=str_wrap(labels[it],width=width), las=1, cex=.25, tick=F, pos=toplot[[j]][[topicid[i]]]$cis[1])
            it = it + 1
          }
        }
      }
    }
}

###################
#Difference method#
###################
if(method=="difference"){
  if(is.null(cov.value1) | is.null(cov.value2)) stop("Cov.value1 and cov.value2 must be specified for categorical variables using method of difference")
  if(length(cov.value1)!=length(cov.value2)) stop("Cov.value1 and cov.value2 must be the same length")
  if(class(cov.value1)=="character"){
    lev <- levels(object$data[,covariate])
  }
  
  #Estimate simulated values to plot
  toplot <- list()
  for(j in 1:length(cov.value1)){
    toplot[[j]] <- list()
    for(i in 1:length(object$topics)){
      if(class(cov.value1)=="character"){
        treat.id <- which(lev==cov.value1[j])
        ct <- contr.treatment(lev)
        cvector[covid] <- ct[rownames(ct)==cov.value1[j]]
        if(any(parse)){
          if(covint>0 & is.null(int.value)){
            cvector[intterms+1] <- ct[rownames(ct)==cov.value1[j]]*median(xmat[,which(matnames==which(labels1==intlabs[cint]))])
          }
          if(covint>0 & !is.null(int.value)){
            cvector[intterms+1] <- ct[rownames(ct)==cov.value1[j]]*int.value
          }
        }
        simst <- simbetas[[i]]%*%cvector
        cvector[covid] <- ct[rownames(ct)==cov.value2[j]]
        if(any(parse)){
          if(covint>0 & is.null(int.value)){
            cvector[intterms+1] <- ct[rownames(ct)==cov.value2[j]]*median(xmat[,which(matnames==which(labels1==intlabs[cint]))])
          }
          if(covint>0 & !is.null(int.value)){
            cvector[intterms+1] <- ct[rownames(ct)==cov.value2[j]]*int.value
          }
        }
        simsc <- simbetas[[i]]%*%cvector
      }
      if(class(cov.value1)!="character"){
        cvector[covid] <- cov.value1[j]
        simst <- simbetas[[i]]%*%cvector
        if(any(parse)){
          if(covint>0 & is.null(int.value)){
            cvector[intterms+1] <- cov.value1[j]*median(xmat[,which(matnames==which(labels1==intlabs[cint]))])
          }
          if(covint>0 & !is.null(int.value)){
            cvector[intterms+1] <- cov.value1[j]*int.value
          }

        }
        cvector[covid] <- cov.value2[j]
        if(any(parse)){
          if(covint>0 & is.null(int.value)){
            cvector[intterms+1] <- cov.value2[j]*median(xmat[,which(matnames==which(labels1==intlabs[cint]))])
          }
          if(covint>0 & !is.null(int.value)){
            cvector[intterms+1] <- cov.value2[j]*int.value
          }
        }
        simsc <- simbetas[[i]]%*%cvector
      }
      toplot[[j]][[i]] <- list(diff=mean(simst-simsc), cis=quantile(simst-simsc, c(offset,1-offset)))
    }
  }
  
  #Adapt labels for difference method
  newlabels <- list()      
  for(i in 1:length(cov.value1)){
    newlabels[[i]] <- vector(length=length(object$topics))
    for(k in 1:length(object$topics)){
      if(verbose.labels){
        newlabels[[i]][k] <- paste(labels[k], " (Covariate Level ", cov.value1[i], " Compared to ", cov.value2[i], ")", sep="")
      } else {
        newlabels[[i]][k] <- paste(labels[k])
      }
    }
  }

  #Plot estimates
  topicid <- rev(which(object$topics%in%topics))
  if (is.null(xlim)) if(labeltype!="numbers") xlim <- c(min(unlist(toplot)) - 1.5*abs(min(unlist(toplot))), max(unlist(toplot)))
  if (is.null(xlim)) if(labeltype=="numbers") xlim <- c(min(unlist(toplot)) - .2*abs(min(unlist(toplot))), max(unlist(toplot)))
  if (is.null(ylim)) ylim <- c(0,length(topicid)*length(cov.value1)+1)
  if (is.null(ylab)) ylab <- ""
  if (add==F) plot(0,0, col="white", xlim=xlim, ylim=ylim, yaxt="n", ylab=ylab, xlab=xlab, main=main, family=family,...)
  lines(c(0,0), c(0, length(topicid)*length(cov.value1)+2), lty=2)
  it = 1
  for(j in 1:length(cov.value1)){
    for(i in 1:length(topicid)){
      points(toplot[[j]][[topicid[i]]]$diff,(j-1)*length(topicid)+i, pch=16)
      lines(c(toplot[[j]][[topicid[i]]]$cis[1], toplot[[j]][[topicid[i]]]$cis[2]), c((j-1)*length(topicid)+i,(j-1)*length(topicid)+i))
      #text(toplot[[j]][[topicid[i]]]$cis[1],(j-1)*length(topicid)+i, str_wrap(newlabels[[j]][topicid[i]],width=width), cex=1, family=family)      
      if(labeltype!="custom") axis(2, at=(j-1)*length(topicid)+i, labels=str_wrap(newlabels[[j]][topicid[i]],width=width), las=1, cex=.25, tick=F, pos=toplot[[j]][[topicid[i]]]$cis[1], family=family)
      if(labeltype=="custom"){
        axis(2, at=(j-1)*length(topicid)+i, labels=str_wrap(labels[it],width=width), las=1, cex=.25, tick=F, pos=toplot[[j]][[topicid[i]]]$cis[1], family=family)
        it = it + 1
      }
    }
  }
}
####################################
#Continuous method, no interactions#
####################################
if((method=="continuous" | method=="spline") & !any(parse)){
  topicid <- which(object$topics%in%topics)
  
  #Estimate simulated values to plot
  toplot <- list()
  if(method=="continuous"){
    uvals <- seq(min(xmat[,covid]), max(xmat[,covid]), length=npoints)
  }
  if(method=="spline"){
    minimum<- min(object$data[,c(labels1[which(term.labels==covspline)])])
    maximum <- max(object$data[,c(labels1[which(term.labels==covspline)])])
    uvals <- seq(minimum, maximum, length=npoints)                   
  }
  
  toplot <- list()
  for(j in 1:length(object$topics)){
    toplot[[j]] <- list()
    toplot[[j]]$mean <- NULL
    toplot[[j]]$cis <- matrix(nrow=0, ncol=2)
    for(i in 1:length(uvals)){
      if(method=="continuous"){
        cvector[covid] <- uvals[i]
      }
      if(method=="spline"){
        cvector[covid] <- predict(object$modelframe[[covspline]],uvals[i]) 
      }
      sims <- simbetas[[j]]%*%cvector
      toplot[[j]]$mean <- c(toplot[[j]]$mean, mean(sims))
      toplot[[j]]$cis <- rbind(toplot[[j]]$cis, quantile(sims, c(offset, 1-offset)))
    }
  }

  #Plot estimates
  if (is.null(xlim)) xlim <- c(min(uvals), max(uvals))
  if (is.null(ylim)) ylim <- c(min(unlist(toplot), na.rm=T), max(unlist(toplot), na.rm=T))
  if (is.null(ylab)) ylab <- "Expected Topic Proportion"
  if (add==F) plot(0, 0,col="white",xlim=xlim, ylim=ylim, main=main, xlab=xlab, ylab=ylab,  ...)
  if (is.null(linecol)) cols = rainbow(length(topicid))
  if (!is.null(linecol)) cols=linecol
  for(i in 1:length(topicid)) {
    lines(uvals, toplot[[topicid[i]]]$mean, col=cols[i])
    lines(uvals, toplot[[topicid[i]]]$cis[,1], col=cols[i],lty=2)
    lines(uvals, toplot[[topicid[i]]]$cis[,2], col=cols[i],lty=2)
  }
  if(printlegend){
    if(labeltype!="custom") legend(xlim[1], ylim[2], labels[topicid], cols)
    if(labeltype=="custom") legend(xlim[1], ylim[2], labels, cols)
  }
}

#################################
#Continuous method, interactions#
#################################
if((method=="continuous" | method=="spline") & any(parse)){
  if(!is.null(int.value)){
    warning("int.value does not affect continuous interactions")
  }
  if(length(unique(xmat[,which(matnames==which(labels1==intlabs[cint]))]))>2){
    stop("For method continuous with interaction, interacted variable must be binary.")
  }
  if(length(topics)>1){
    stop("Interaction plot can only be used for one topic at a time.")
  }

  interactlevels <- unique(xmat[,which(matnames==which(labels1==intlabs[cint]))])
  
  cvector[which(matnames==which(labels1==intlabs[cint]))] <- interactlevels[1]
  cvector2 <- cvector
  cvector2[which(matnames==which(labels1==intlabs[cint]))] <- interactlevels[2]
  
  topicid <- which(object$topics%in%topics)
  #Estimate simulated values to plot
  toplot <- list()
  toplot2 <- list()
  
  if(method=="continuous"){
    uvals <- seq(min(xmat[,covid]), max(xmat[,covid]), length=npoints)
  }
  if(method=="spline"){
    stop("Interaction plots not yet supported for splines")
    minimum<- min(object$data[,c(labels1[which(term.labels==covspline)])])
    maximum <- max(object$data[,c(labels1[which(term.labels==covspline)])])
    uvals <- seq(minimum, maximum, length=npoints)                   
  }

  for(j in 1:length(object$topics)){
    toplot[[j]] <- list()
    toplot[[j]]$mean <- NULL
    toplot[[j]]$cis <- matrix(nrow=0, ncol=2)
    toplot2[[j]] <- list()
    toplot2[[j]]$mean <- NULL
    toplot2[[j]]$cis <- matrix(nrow=0, ncol=2)

    for(i in 1:length(uvals)){
      if(method=="continuous"){
        cvector[covid] <- cvector2[covid] <- uvals[i]
        cvector[intterms+1] <- uvals[i]*interactlevels[1]
        cvector2[intterms+1] <- uvals[i]*interactlevels[2]
      }
      if(method=="spline"){
        cvector[covid] <- cvector2[covid] <- predict(object$modelframe[[covspline]],uvals[i])
        cvector[which(matnames==which(intterms==T))] <- predict(object$modelframe[[covspline]],uvals[i])*interactlevels[1]
        cvector[which(matnames==which(intterms==T))] <- predict(object$modelframe[[covspline]],uvals[i])*interactlevels[2]
      }
      sims <- simbetas[[j]]%*%cvector
      sims2 <- simbetas[[j]]%*%cvector2
      toplot[[j]]$mean <- c(toplot[[j]]$mean, mean(sims))
      toplot[[j]]$cis <- rbind(toplot[[j]]$cis, quantile(sims, c(offset, 1-offset)))
      toplot2[[j]]$mean <- c(toplot2[[j]]$mean, mean(sims2))
      toplot2[[j]]$cis <- rbind(toplot2[[j]]$cis, quantile(sims2, c(offset, 1-offset)))
      
    }
  }

  #Plot estimates
  if (is.null(xlim)) xlim <- c(min(uvals), max(uvals))
  if (is.null(ylim)) ylim <- c(min(c(min(unlist(toplot), na.rm=T),min(unlist(toplot2), na.rm=T))), max(c(max(unlist(toplot), na.rm=T),max(unlist(toplot2), na.rm=T))))
  if (is.null(ylab)) ylab <- "Expected Topic Proportion"
  if (add==F) plot(0, 0,col="white",xlim=xlim, ylim=ylim, main=main, xlab=xlab, ylab=ylab,  ...)
  if (is.null(linecol)) cols = rainbow(2)
  if (!is.null(linecol)) cols=linecol
  for(i in 1:length(topicid)) {
    lines(uvals, toplot[[topicid[i]]]$mean, col=cols[1])
    lines(uvals, toplot[[topicid[i]]]$cis[,1], col=cols[1],lty=2)
    lines(uvals, toplot[[topicid[i]]]$cis[,2], col=cols[1],lty=2)
    lines(uvals, toplot2[[topicid[i]]]$mean, col=cols[2])
    lines(uvals, toplot2[[topicid[i]]]$cis[,1], col=cols[2],lty=2)
    lines(uvals, toplot2[[topicid[i]]]$cis[,2], col=cols[2],lty=2)
  }
  if(printlegend){
    if(labeltype!="custom")   legend(xlim[1], ylim[2], c(paste(labels[topicid], intlabs[cint], "=", interactlevels[1]),paste(labels[topicid], intlabs[cint], "=", interactlevels[2])) , cols)
    if(labeltype=="custom") legend(xlim[1], ylim[2], labels, cols)
  }
}
}


#This function uses a categorical covariate and plots the means of the thetas. 
#It should eventually be subsummed.
ploteffect.catmeans<-function(model,covariate=covariate,sims=500,k=k,xlab=xlab,ylab=ylab,main=main) {
theta <- thetapost.global(model, sims)

k <- k
out <- list()
rdraw <- function(vec) {
  x <- mean(vec)
  sig <- sd(vec)/sqrt(length(vec))
  rnorm(100,x,sig)
}
z<-max(covariate)
y<-min(covariate)
for(g in y:z) {
  index <- which(covariate==g)
  ndocs <- length(index)
  group <- theta[index]
  rm(index)
  ndocs <- length(group)
  group <- do.call(rbind, group)
  index <- rep(1:sims, times=ndocs)
  group <- split(group, index)
  group <- lapply(group, matrix, nrow=ndocs)
  draws <- unlist(lapply(group, function(x) rdraw(x[,k])))
  out[[g+1]] <- draws
}

xmeans <- unlist(lapply(out,mean))
xlower <- unlist(lapply(out, quantile, .025))
xupper <- unlist(lapply(out, quantile, .975))

plot(0:z, xmeans,ylim=c(min(xlower), max(xupper)), xlab=xlab, ylab=ylab, main=paste(c("Topic",k),collapse=" "))
segments(0:z, xlower, 0:z, xupper,lty=2)
}


