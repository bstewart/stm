#' Exclusivity
#' 
#' Calculate an exclusivity metric for an STM model.
#' 
#' In Roberts et al 2014 we proposed using the Mimno et al 2011 \code{\link{semanticCoherence}} metric 
#' for helping with topic model selection. We found that semantic coherence alone is relatively easy to
#' achieve by having only a couple of topics which all are dominated by the most common words.  Thus we
#' also proposed an exclusivity measure.  
#' 
#' Our exclusivity measure includes some information on word frequency as well.  It is based on the FREX
#' labeling metric (\code{\link{calcfrex}}) with the weight set to .7 in favor of exclusivity by default.
#'  
#'  This function is currently marked with the keyword internal because it does not have much error checking.
#'
#' @param model the STM object
#' @param M the number of top words to consider per topic
#' @param frexw the frex weight
#' 
#' @return a numeric vector containing exclusivity for each topic
#' 
#' @references 
#' Mimno, D., Wallach, H. M., Talley, E., Leenders, M., & McCallum, A. (2011, July). 
#' "Optimizing semantic coherence in topic models." In Proceedings of the Conference on Empirical Methods in 
#' Natural Language Processing (pp. 262-272). Association for Computational Linguistics. Chicago
#' 
#' Bischof and Airoldi (2012) "Summarizing topical content with word frequency and exclusivity"
#' In Proceedings of the International Conference on Machine Learning.
#' 
#' Roberts, M., Stewart, B., Tingley, D., Lucas, C., Leder-Luis, J., Gadarian, S., Albertson, B., et al. (2014). 
#' "Structural topic models for open ended survey responses." American Journal of Political Science, 58(4), 1064-1082.
#' @seealso \code{\link{searchK}} \code{\link{plot.searchK}} \code{\link{semanticCoherence}}
#' @keywords internal
#' @examples 
#' exclusivity(gadarianFit)
#' @export
exclusivity <- function(model, M=10, frexw=.7){
  w <- frexw
  if(length(model$beta$logbeta)!=1) stop("Exclusivity calculation only designed for models without content covariates")
  tbeta <- t(exp(model$beta$logbeta[[1]]))
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