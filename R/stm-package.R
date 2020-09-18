#' Structural Topic Model
#' 
#' This package implements the Structural Topic Model, a general approach to
#' including document-level metadata within mixed-membership topic models. To
#' read the vignette use \code{vignette('stmVignette')}.
#' 
#' Functions to manipulate documents: \code{\link{textProcessor}}
#' \code{\link{readCorpus}} \code{\link{prepDocuments}}
#' 
#' Functions to fit the model: \code{\link{stm}} \code{\link{selectModel}}
#' \code{\link{manyTopics}} \code{\link{searchK}}
#' 
#' Functions to summarize a model: \code{\link{labelTopics}}
#' \code{\link{summary.STM}} \code{\link{findThoughts}}
#' 
#' Functions for Post-Estimation: \code{\link{estimateEffect}}
#' \code{\link{topicCorr}} \code{\link{permutationTest}}
#' 
#' Plotting Functions: \code{\link{plot.STM}} \code{\link{plot.estimateEffect}}
#' \code{\link{plot.topicCorr}} \code{\link{plot.STMpermute}}
#' \code{\link{plotQuote}} \code{\link{plotTopicLoess}}
#' \code{\link{plotModels}} \code{\link{topicQuality}}
#' 
#' Pre-Fit Models and Data: \code{\link{gadarian}} \code{\link{gadarianFit}}
#' \code{\link{poliblog5k}}
#' 
#' 
#' @name stm-package
#' @docType package
#' @author Author: Margaret E. Roberts, Brandon M. Stewart and Dustin Tingley
#' 
#' Maintainer: Brandon Stewart <bms4@@princeton.edu>
#' @seealso \code{\link{stm}}
#' @references Roberts, M., Stewart, B., Tingley, D., and Airoldi, E. (2013)
#' "The structural topic model and applied social science." In Advances in
#' Neural Information Processing Systems Workshop on Topic Models: Computation,
#' Application, and Evaluation.
#' 
#' Roberts, M., Stewart, B., Tingley, D., Lucas, C., Leder-Luis, J., Gadarian,
#' S., Albertson, B., Albertson, B. and Rand, D. (2014). "Structural topic
#' models for open ended survey responses." American Journal of Political
#' Science.
#' 
#' Additional papers at: structuraltopicmodel.com
#' @keywords package
#' 
#' @import Matrix
#' @importFrom Rcpp evalCpp
#' @importFrom graphics abline axis hist legend lines par plot points segments smoothScatter text title
#' @importFrom stats aggregate as.formula coef cor cov lm loess median model.frame model.response model.matrix na.omit optim optimize pchisq predict quantile rbinom rgamma rnorm runif terms
#' @useDynLib stm, .registration = TRUE
NULL

#' Gadarian and Albertson data
#' 
#' This data set
#' contains variables from Gadarian and Albertson (2014). The experiment had
#' those in the treatment condition write about what made them anxious about
#' immigration. The control condition just had subjects write about
#' immigration.
#' 
#' 
#' @name gadarian
#' @aliases gadarian gadarianFit
#' @docType data
#' @format A data frame with 351 observations on the following 3 variables.
#' \describe{ 
#' \item{\code{MetaID}}{A numeric vector containing identification
#' numbers; not used for analysis} 
#' \item{\code{treatment}}{A numeric vector
#' indicating treatment condition} 
#' \item{\code{pid_rep}}{A numeric vector of
#' party identification} 
#' \item{\code{open.ended.response}}{A character vector
#' of the subject's open ended response} 
#' }
#' @source Gadarian, Shana Kushner, and Bethany Albertson. "Anxiety,
#' immigration, and the search for information." Political Psychology 35.2
#' (2014): 133-164.
#' 
#' Roberts, Margaret E., Brandon M. Stewart, Dustin Tingley, Christopher Lucas,
#' Jetson Leder-Luis, Shana Kushner Gadarian, Bethany Albertson, and David G.
#' Rand.  "Structural Topic Models for Open-Ended Survey Responses." American
#' Journal of Political Science 58, no 4 (2014): 1064-1082.
#' @keywords datasets
#' @examples
#' \donttest{
#' 
#' head(gadarian)
#' #Process the data for analysis.
#' temp<-textProcessor(documents=gadarian$open.ended.response,metadata=gadarian)
#' meta<-temp$meta
#' vocab<-temp$vocab
#' docs<-temp$documents
#' out <- prepDocuments(docs, vocab, meta)
#' docs<-out$documents
#' vocab<-out$vocab
#' meta <-out$meta
#' }
NULL

#' CMU 2008 Political Blog Corpus
#' 
#' A 5000 document sample from CMU 2008 Political Blog Corpus (Eisenstein and
#' Xing 2010).  Blog posts from 6 blogs during the U.S. 2008 Presidential
#' Election.
#' 
#' This is a random sample of the larger CMU 2008 Political Blog Corpus
#' collected by Jacob Eisenstein and Eric Xing.  Quoting from their
#' documentation: "[The blogs] were selected by the following criteria: the
#' Technorati rankings of blog authority, ideological balance, coverage for the
#' full year 2008, and ease of access to blog archives. In the general election
#' for U.S. President in 2008, the following blogs supported Barack Obama:
#' Digby, ThinkProgress, and Talking Points Memo. John McCain was supported by
#' American Thinker, Hot Air, and Michelle Malkin. In general, the blogs that
#' supported Obama in the election tend to advocate for similar policies and
#' candidates as the Democratic party; and the blogs that supported McCain tend
#' to advocate Republican policies and candidates. Digby, Hot Air and Michelle
#' Malkin are single-author blogs; the others have multiple authors."
#' 
#' @name poliblog5k
#' @aliases poliblog5k poliblog5k.docs poliblog5k.voc poliblog5k.meta
#' @docType data
#' @format A data frame with 5000 observations on the following 4 variables.
#' \describe{ 
#' \item{\code{rating}}{a factor variable giving the partisan
#' affiliation of the blog (based on who they supported for president)}
#' \item{\code{day}}{the day of the year (1 to 365).  All entries are from
#' 2008.} 
#' \item{\code{blog}}{a two digit character code corresponding to the
#' name of the blog. They are: American Thinker (at), Digby (db), Hot Air (ha),
#' Michelle Malkin (mm), Think Progress (tp), Talking Points Memo (tpm)}
#' \item{\code{text}}{the first 50 characters (rounded to the nearest full
#' word).} 
#' }
#' @source Jacob Eisenstein and Eric Xing (2010) "The CMU 2008 Political Blog
#' Corpus." Technical Report Carnegie Mellon University.
#' http://sailing.cs.cmu.edu/socialmedia/blog2008.html
#' @keywords datasets
#' @examples
#' 
#' \donttest{
#' 
#' data(poliblog5k)
#' head(poliblog5k.meta)
#' head(poliblog5k.voc)
#' 
#' stm1 <- stm(poliblog5k.docs, poliblog5k.voc, 3,
#' prevalence=~rating, data=poliblog5k.meta)
#' }
#' 
NULL