#Outputs median or mode of a dataset depending on the type
#type: class of variable
#data: data in which the variable exists
outputdata <- function(type, data){
  if(type=="character"|type=="logical"){
    tab <- table(data[[names(type)]])
    mode <- names(tab)[which.max(tab)]
    return(factor(mode, levels=levels(as.factor(data[[names(type)]]))))
  }
  if(type=="factor"){
    tab <- table(data[,names(type)])
    mode <- names(tab)[which.max(tab)]
    return(factor(mode, levels=levels(data[[names(type)]])))
  }
  if(type=="numeric"|type=="integer"){
    return(median(data[[names(type)]]))
  }
}

#Takes method for estimateEffect and produces a control matrix to
#simulate from
#prep: output of estimateEffect
#covariate: covariate of interest
#method: pointestimate, difference, or continuous
#cov.value1 and cov.value2: for method=difference
#npoints: for the continuous covariate, how many points from min to
#max to draw from.
produce_cmatrix <- function(prep, covariate, method,cov.value1=NULL,
                            cov.value2=NULL, npoints=100, moderator=NULL, moderator.value=NULL){

  #Find type of each variable
  types <- lapply(prep$data, function(x) class(x)[1])
  types <- unlist(types[prep$varlist]) 
  #switched the below out.
  #types <- sapply(prep$varlist,function (x) class(prep$data[,x]))

  #Make control matrix
  #What is the covariate of interest and what are the controls?
  covariateofinterest <- which(prep$varlist==covariate)
  controls <- which(prep$varlist!=covariate)
  
  if(method=="pointestimate"){    
                                        #Start cdata with the variable of interest
    if(types[covariateofinterest]=="character") cdata <- data.frame(factor(unique(prep$data[[covariate]])))
      
    if(types[covariateofinterest]=="factor") cdata <- data.frame(unique(prep$data[[covariate]]))
    
    if(types[covariateofinterest]=="numeric" |
       types[covariateofinterest]=="integer") cdata <-
       data.frame(unique(prep$data[[covariate]]))
    
    names(cdata) <- covariate
  }
  if(method=="difference"){
    if(types[covariateofinterest]=="character" |
       types[covariateofinterest]=="factor") {
      lev <- levels(as.factor(prep$data[[covariate]]))
      x <- c(as.character(cov.value1), as.character(cov.value2))
      cdata <- base::data.frame(factor(x,levels=lev))
      colnames(cdata) <- covariate 
      rm(x,lev)
    }
    if(types[covariateofinterest]=="numeric" |
       types[covariateofinterest]=="integer") cdata <-  base::data.frame(c(cov.value1, cov.value2))
    names(cdata) <- covariate
  }
  
  if(method=="continuous"){
    if(types[covariateofinterest]=="character" |
       types[covariateofinterest]=="factor")
    stop("Covariate of interest must be numeric")
    if(types[covariateofinterest]=="numeric" |
       types[covariateofinterest]=="integer") cdata <-
  base::data.frame(seq(min(prep$data[[covariate]]),max(prep$data[[covariate]]),
    length.out=npoints))
    names(cdata) <- covariate
  }
                                        #Insert the values of the controls  
  if(length(controls)>0){
      for(i in 1:length(controls)){
        cdata[,prep$varlist[controls[i]]] <-outputdata(types[controls[i]], prep$data)
      }
    }

  #Insert the value for the interaction, if applicable
  if(!is.null(moderator) & !is.null(moderator.value)){
    if(is.factor(prep$data[[moderator]])|is.character(prep$data[[moderator]])){
      cdata[,moderator] <- factor(moderator.value, levels=levels(as.factor(prep$data[[moderator]])))
    }else{
      cdata[,moderator] <- moderator.value
    }
  }
  if(!is.null(moderator) & is.null(moderator.value)){
    stop("Please specify the value of the moderator")
  }
  if(is.null(moderator) & !is.null(moderator.value)){
    stop("Please specify the moderator")
  }
  #Reorder to reflect original data
  if(ncol(cdata)>1){
    cdata <- cdata[,names(prep$data)]
  }

  #Get model.matrix
  #cmatrix <- parseFormulas(prep, cdata)
  cmatrix <- makeDesignMatrix(prep$formula, prep$data, cdata, sparse=FALSE)
  return(list(cdata=cdata,cmatrix=cmatrix))
}
