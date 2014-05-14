#function to choose the pareto dominant model, or randomly choose among model candidates that are not weakly dominanted by other models
#utility function for manyTopics

paretosingle <- function(z) {
  
  m<-matrix(NA,nrow=length(z$semcoh),ncol=2)
  
  for(i in 1:nrow(m)){
    
    m[i,1]<-as.numeric(mean(unlist(z$semcoh[i])))
    
    if(!z$exclusivity[[1]][1]=="Exclusivity not calculated for models with content covariates"){
      m[i,2]<-as.numeric(mean(unlist(z$exclusivity[i])))
    } else {
      stop("manyTopics function not yet designed for models with content variable.") 
    }
  }
 
  d1max <- max(m[,1])
  d2max <- max(m[,2])
  weakcandidates <- m[,1]==d1max | m[,2]==d2max
  strongcandidates <- m[,1]==d1max & m[,2]==d2max
  s = which(strongcandidates)
  w = which(weakcandidates)
  if (length(s)>0) {
    x = s
  }
  else {
    x = w
  }
  if (length(x)==1) {return(x)}
  else {
    return(sample(x,size=1))
  }
}


#loop over many different numbers of topics.
manyTopics <- function(documents, vocab, K, prevalence, content, 
                       data = NULL,max.em.its = 100, verbose = TRUE, 
                       init.type = "LDA", 
                       emtol = 1e-05, seed = NULL, runs = 50, 
                       frexw = 0.7, net.max.em.its = 2, 
                       netverbose = FALSE, M = 10,...) {
  out<-list()
  semcoh<-exclusivity<-list()
  for(i in 1:length(K)) {
    
    models<- selectModel(documents, vocab, K[i], prevalence, content, data = data, 
                         max.em.its = max.em.its, verbose = verbose, init.type = init.type, emtol = emtol, seed = seed, runs = runs, 
                         frexw = frexw, net.max.em.its = net.max.em.its, netverbose = netverbose, M = M, 
                         ...)
    j<-paretosingle(models)

    out[[i]]<-models$runout[[j]]
    exclusivity[[i]]<-models$exclusivity[[j]]
    semcoh[[i]]<-models$semcoh[[j]]
    j<-NULL
  }
  
  toreturn<-list(out=out,exclusivity=exclusivity,semcoh=semcoh)
  return(toreturn)
}




