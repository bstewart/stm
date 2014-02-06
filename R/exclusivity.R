#rewritten exclusivity code 1/25/2014 to make it a bit faster and add some errors
# benchmarked and got exact same answers against old code.

exclusivity <- function(mod.out, M=10, frexw=.7){
  w <- frexw
  if(length(mod.out$beta$logbeta)!=1) stop("Exclusivity calculation only designed for models without content covariates")
  tbeta <- t(exp(mod.out$beta$logbeta[[1]]))
  s <- rowSums(tbeta)
  mat <- tbeta/s #normed by columns of beta now.

  ex <- apply(mat,2,rank)/nrow(mat)
  fr <- apply(tbeta,2,rank)/nrow(mat)
  frex<- 1/(w/ex + (1-w)/fr)
  index <- apply(tbeta, 2, order, decreasing=TRUE)[1:M,]
  out <- vector(length=ncol(tbeta)) 
  for(i in 1:ncol(frex)) {
    out[i] <- sum(frex[index[,i],i])
  }
  out
}