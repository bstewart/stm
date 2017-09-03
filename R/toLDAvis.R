#' Wrapper to launch LDAvis topic browser.
#' 
#' Tool for exploring topic/word distributions using LDAvis topic browser.
#' 
#' Tool for exploring topic/word distributions using LDAvis topic browser.
#' Development build of LDAvis available at https://github.com/cpsievert/LDAvis
#' or download from CRAN. Note: LDAvis may renumber the topics. 
#' 
#' @param mod STM object. Output from stm function.
#' @param docs Documents object passed to \code{stm} in this package's standard
#' format (see the documentation in \code{\link{stm}}.
#' @param R Passed to \code{\link[LDAvis]{createJSON}} "integer, the number of
#' terms to display in the barcharts of the interactive viz. Default is 30.
#' Recommended to be roughly between 10 and 50."
#' @param plot.opts Passed to \code{\link[LDAvis]{createJSON}} "a named list
#' used to customize various plot elements. By default, the x and y axes are
#' labeled 'PC1' and 'PC2' (principal components 1 and 2), since jsPCA is the
#' default scaling method. "
#' @param lambda.step Passed to \code{\link[LDAvis]{createJSON}} "a value
#' between 0 and 1. Determines the interstep distance in the grid of lambda
#' values over which to iterate when computing relevance. Default is 0.01.
#' Recommended to be between 0.01 and 0.1."
#' @param out.dir Passed to \code{\link[LDAvis]{serVis}} "directory to store
#' html/js/json files."
#' @param open.browser Passed to \code{\link[LDAvis]{serVis}} "Should R open a
#' browser? If yes, this function will attempt to create a local file server
#' via the servr package. This is necessary since the javascript needs to
#' access local files and most browsers will not allow this."
#' @param as.gist Passed to \code{\link[LDAvis]{serVis}} "should the vis be
#' uploaded as a gist? Will prompt for an interactive login if the GITHUB_PAT
#' environment variable is not set"
#' @references Carson Sievert and Kenny Shirley. LDAvis: Interactive
#' Visualization of Topic Models. R package version 0.3.1.
#' https://github.com/cpsievert/LDAvis
#' @examples
#' 
#' \dontrun{
#' 
#' mod <- stm(poliblog5k.docs, poliblog5k.voc, K=25,
#'            prevalence=~rating, data=poliblog5k.meta,
#'            max.em.its=2, init.type="Spectral") 
#' #please don't run a model with 2 iterations
#' #this is done here to make it run quickly.
#' toLDAvis(mod=mod, docs=poliblog5k.docs)
#' 
#' }
#' @export
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
LDAvis::serVis(f,out.dir=out.dir,open.browser=open.browser,as.gist=as.gist,R = R)
}
