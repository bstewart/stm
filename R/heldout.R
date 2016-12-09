
make.heldout <- function(documents, vocab, N=floor(.1*length(documents)), 
                         proportion=.5, seed=NULL) {
  if(!is.null(seed)) set.seed(seed)
  index <- sort(sample(1:length(documents), N))
  pie <- proportion
  missing <- vector(mode="list", length=N)
  ct <- 0
  for(i in index) {
    ct <- ct + 1
    doc <- documents[[i]]  
    if(ncol(doc)<2) next
    doc <- rep(doc[1,], doc[2,])
    #how many tokens to sample? The max ensures at least one is sampled
    nsamp <- max(1,floor(pie*length(doc)))
    ho.index <- sample(1:length(doc), nsamp)
    tab <- tabulate(doc[ho.index])
    missing[[ct]] <- rbind(which(tab>0), tab[tab>0])
    tab <- tabulate(doc[-ho.index])
    documents[[i]] <- rbind(which(tab>0), tab[tab>0])
  }
  missing <- list(index=index, docs=missing)
  #check the vocab
  indices <- sort(unique(unlist(lapply(documents, "[", 1, ))))
  
  #all sorts of nonsense ensues if there is missingness
  if(length(indices)!=length(vocab)) {
    remove <- which(!(1:length(vocab)%in% indices))
    newind <- rep(0, length(vocab))
    newind[indices] <- 1:length(indices)
    new.map <- cbind(1:length(vocab), newind)
    #renumber the missing elements and remove 0's
    missing$docs <- lapply(missing$docs, function(d) {
      d[1,] <- new.map[match(d[1,], new.map[,1]),2]
      return(d[,d[1,]!=0, drop=FALSE])
    })
    #apply the same process to the documents
    documents <- lapply(documents, function(d) {
      d[1,] <- new.map[match(d[1,], new.map[,1]),2]
      return(d[,d[1,]!=0, drop=FALSE])
    })
    
    lens <- unlist(lapply(missing$docs, length))
    if(any(lens==0)) {
      missing$docs <- missing$docs[lens!=0]
      missing$index <- missing$index[lens!=0]
    }
    vocab <- vocab[indices]
  }
  #hooray.  return some stuff.
  return(list(documents=documents,vocab=vocab, missing=missing))
}

eval.heldout <- function(model, missing) {
  heldout <- vector(length=length(missing$index))
  ntokens <- vector(length=length(missing$index))
  beta <- lapply(model$beta$logbeta, exp)
  bindex <- model$settings$covariates$betaindex[missing$index]
  for(i in 1:length(missing$index)) {
    docid <- missing$index[i]
    words <- missing$docs[[i]][1,]
    probs <- model$theta[docid,]%*%beta[[bindex[i]]][,words]
    probs <- rep(probs, missing$docs[[i]][2,])
    heldout[i] <- mean(log(probs)) 
    ntokens[i] <- sum(words)
  }
  out <- list(expected.heldout=mean(heldout, na.rm=TRUE), doc.heldout=heldout,
              index=missing$index, ntokens=ntokens) #the mean here to deal with 0 length docs
  return(out)
}




