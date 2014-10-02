semanticCoherence <- function(model.out, documents, M){
  if(length(model.out$beta$logbeta)!=1) {
    result <- 0
    for(i in 1:length(model.out$beta$logbeta)){
      subset <- which(model.out$settings$covariates$betaindex==i)
      triplet <- doc.to.ijv(documents[subset])
      mat <- simple_triplet_matrix(triplet$i, triplet$j,triplet$v, ncol=model.out$settings$dim$V)
      result = result + semCoh1beta(mat, M, beta=model.out$beta$logbeta[[i]])*length(subset)
    }
    return(result/length(documents))
  }
  else {
    beta <- model.out$beta$logbeta[[1]]
    #Get the Top N Words
    top.words <- apply(beta, 1, order, decreasing=TRUE)[1:M,]
    wordlist <- unique(as.vector(top.words))
    triplet <- doc.to.ijv(documents)
    mat <- simple_triplet_matrix(triplet$i, triplet$j,triplet$v, ncol=model.out$settings$dim$V)
    result = semCoh1beta(mat, M, beta=beta)
  return(result)
  }
}


semCoh1beta <- function(mat, M, beta){
   #Get the Top N Words
  top.words <- apply(beta, 1, order, decreasing=TRUE)[1:M,]
  wordlist <- unique(as.vector(top.words))
  mat <- mat[,wordlist]
  mat$v <- ifelse(mat$v>1, 1,mat$v) #binarize
  
  #do the cross product to get co-occurences
  cross <- tcrossprod_simple_triplet_matrix(t(mat))
  
  #create a list object with the renumbered words (so now it corresponds to the rows in the table)
  temp <- match(as.vector(top.words),wordlist)
  labels <- split(temp, rep(1:nrow(beta), each=M))
  
  #Note this could be done with recursion in an elegant way, but let's just be simpler about it.
  sem <- function(ml,cross) {
    m <- ml[1]; l <- ml[2]
    log(.01 + cross[m,l]) - log(cross[l,l] + .01)
  }
  result <- vector(length=nrow(beta))
  for(k in 1:nrow(beta)) {
    grid <- expand.grid(labels[[k]],labels[[k]])
    colnames(grid) <- c("m", "l") #corresponds to original paper
    grid <- grid[grid$m > grid$l,]
    calc <- apply(grid,1,sem,cross)
    result[k] <- sum(calc)
  }
  return(result)
}
