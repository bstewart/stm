#' Semantic Coherence
#' 
#' Calculate semantic coherence (Mimno et al 2011) for an STM model.
#' 
#'  Semantic coherence is a metric related to pointwise mutual information that was introduced
#'  in a paper by David Mimno, Hanna Wallach and colleagues (see references),  The paper details a series
#'  of manual evaluations which show that their metric is a reasonable surrogate for human judgment.
#'  The core idea here is that in models which are semantically coherent the words which are most
#'  probable under a topic should co-occur within the same document.
#'  
#'  One of our observations in Roberts et al 2014 was that semantic coherence alone is relatively easy to
#'  achieve by having only a couple of topics which all are dominated by the most common words.  Thus we
#'  suggest that users should also consider \code{\link{exclusivity}} which provides a natural counterpoint.
#'  
#'  This function is currently marked with the keyword internal because it does not have much error checking.
#'
#' @param model the STM object
#' @param documents the STM formatted documents
#' @param M the number of top words to consider per topic
#' 
#' @return a numeric vector containing semantic coherence for each topic
#' 
#' @references 
#' Mimno, D., Wallach, H. M., Talley, E., Leenders, M., & McCallum, A. (2011, July). 
#' "Optimizing semantic coherence in topic models." In Proceedings of the Conference on Empirical Methods in 
#' Natural Language Processing (pp. 262-272). Association for Computational Linguistics. Chicago
#' 
#' Roberts, M., Stewart, B., Tingley, D., Lucas, C., Leder-Luis, J., Gadarian, S., Albertson, B., et al. (2014). 
#' "Structural topic models for open ended survey responses." American Journal of Political Science, 58(4), 1064-1082.
#'  http://goo.gl/0x0tHJ		
#' @seealso \code{\link{searchK}} \code{\link{plot.searchK}} \code{\link{exclusivity}}
#' @keywords internal
#' @examples 
#' temp<-textProcessor(documents=gadarian$open.ended.response,metadata=gadarian)
#' meta<-temp$meta
#' vocab<-temp$vocab
#' docs<-temp$documents
#' out <- prepDocuments(docs, vocab, meta)
#' docs<-out$documents
#' vocab<-out$vocab
#' meta <-out$meta
#' set.seed(02138)
#' #maximum EM iterations set very low so example will run quickly.
#' #Run your models to convergence!
#' mod.out <- stm(docs, vocab, 3, prevalence=~treatment + s(pid_rep), data=meta,
#'                max.em.its=5)
#' semanticCoherence(mod.out, docs)
#' @export
semanticCoherence <- function(model, documents, M=10){
  if(!inherits(model, "STM")) stop("model must be an STM object")
                                   
  if(length(model$beta$logbeta)!=1) {
    result <- 0
    for(i in 1:length(model$beta$logbeta)){
      subset <- which(model$settings$covariates$betaindex==i)
      triplet <- doc.to.ijv(documents[subset])
      mat <- slam::simple_triplet_matrix(triplet$i, triplet$j,triplet$v, ncol=model$settings$dim$V)
      result = result + semCoh1beta(mat, M, beta=model$beta$logbeta[[i]])*length(subset)
    }
    return(result/length(documents))
  }
  else {
    beta <- model$beta$logbeta[[1]]
    #Get the Top N Words
    top.words <- apply(beta, 1, order, decreasing=TRUE)[1:M,]
    triplet <- doc.to.ijv(documents)
    mat <- slam::simple_triplet_matrix(triplet$i, triplet$j,triplet$v, ncol=model$settings$dim$V)
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
  cross <- slam::tcrossprod_simple_triplet_matrix(t(mat))
  
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
