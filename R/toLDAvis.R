toLDAvis<-function(mod,docs,R=30,plot.opts=list(xlab ="PC1", ylab = "PC2"),
                    lambda.step=.1,out.dir=tempfile(),open.browser=interactive(),as.gist=FALSE){
if(!requireNamespace("LDAvis",quietly=TRUE)) stop("Please install LDAvis package to use this function. You will also need servr.")
theta<-mod$theta
if(length(mod$beta$logbeta)>1) stop("This function does not yet allow content covariates.")
phi <- exp(mod$beta$logbeta[[1]])
  if(any(phi==0)){
    phi<-phi + .Machine$double.eps
    phi<-phi/rowSums(phi)
  }
vocab <- mod$vocab
doc.length <- as.integer(unlist(lapply(docs, function(x) sum(x[2,]))))
term.frequency <- mod$settings$dim$wcounts$x
f<-LDAvis::createJSON(phi=phi,theta=theta,doc.length=doc.length,vocab=vocab,term.frequency=term.frequency,lambda.step = lambda.step)
LDAvis::serVis(f,out.dir=out.dir,open.browser=open.browser,as.gist=as.gist)
}
