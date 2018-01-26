#' Variational EM for the Structural Topic Model
#' 
#' 
#' Estimation of the Structural Topic Model using semi-collapsed variational
#' EM.  The function takes sparse representation of a document-term matrix, an integer
#' number of topics, and covariates and returns fitted model parameters.
#' Covariates can be used in the prior for topic \code{prevalence}, in the
#' prior for topical \code{content} or both.  See an overview of functions in
#' the package here: \code{\link{stm-package}}
#' 
#' This is the main function for estimating a Structural Topic Model (STM).
#' STM is an admixture with covariates in both mixture components.  Users
#' provide a corpus of documents and a number of topics.  Each word in a
#' document comes from exactly one topic and each document is represented by
#' the proportion of its words that come from each of the K topics.  These
#' proportions are found in the N (number of documents) by K (user specified
#' number of topics) theta matrix.  Each of the K topics are represented as
#' distributions over words.  The K-by-V (number of words in the vocabulary)
#' matrix logbeta contains the natural log of the probability of seeing each
#' word conditional on the topic.
#' 
#' The most important user input in parametric topic models is the number of
#' topics.  There is no right answer to the appropriate number of topics.  More
#' topics will give more fine-grained representations of the data at the
#' potential cost of being less precisely estimated.  The number must be at
#' least 2 which is equivalent to a unidimensional scaling model.  For short
#' corpora focused on very specific subject matter (such as survey experiments)
#' 3-10 topics is a useful starting range.  For small corpora (a few hundred to
#' a few thousand) 5-50 topics is a good place to start.  Beyond these rough
#' guidelines it is application specific.  Previous applications in political
#' science with medium sized corpora (10k to 100k documents) have found 60-100
#' topics to work well.  For larger corpora 100 topics is a useful default
#' size.  Of course, your mileage may vary.
#' 
#' When \code{init.type="Spectral"} and \code{K=0} the number of topics is set
#' using the algorithm in Lee and Mimno (2014).  See vignette for details.  We
#' emphasize here as we do there that this does not estimate the "true" number
#' of topics and does not necessarily have any particular statistical
#' properties for consistently estimating the number of topics.  It can however
#' provide a useful starting point.
#' 
#' The model for topical prevalence includes covariates which the analyst
#' believes may influence the frequency with which a topic is discussed.  This
#' is specified as a formula which can contain smooth terms using splines or by
#' using the function \code{\link{s}}.  The response portion of the formula
#' should be left blank.  See the examples.  These variables can include
#' numeric and factor variables.  While including variables of class
#' \code{Dates} or other non-numeric, non-factor types will work in \code{stm}
#' it may not always work for downstream functions such as
#' \code{\link{estimateEffect}}.
#' 
#' The topical convent covariates are those which affect the way in which a
#' topic is discussed. As currently implemented this must be a single variable
#' which defines a discrete partition of the dataset (each document is in one
#' and only one group).  We may relax this in the future.  While including more
#' covariates in topical prevalence will rarely affect the speed of the model,
#' including additional levels of the content covariates can make the model
#' much slower to converge.  This is due to the model operating in the much
#' higher dimensional space of words in dictionary (which tend to be in the
#' thousands) as opposed to topics.
#' 
#' In addition to the default priors for prevalence, we also make use of the
#' \code{glmnet} package to allow for penalties between the L1 and L2 norm.  In
#' these settings we estimate a regularization path and then select the optimal
#' shrinkage parameter using a user-tuneable information criterion.  By default
#' selecting the \code{L1} option will apply the L1 penalty selecting the
#' optimal shrinkage parameter using AIC. The defaults have been specifically
#' tuned for the STM but almost all the relevant arguments can be changed
#' through the control structure below.  Changing the \code{gamma.enet}
#' parameters allow the user to choose a mix between the L1 and L2 norms.  When
#' set to 1 (as by default) this is the lasso penalty, when set to 0 its the
#' ridge penalty.  Any value in between is a mixture called the elastic net.
#' 
#' The default prior choice for content covariates is now the \code{L1} option.
#' This uses an approximation framework developed in Taddy (2013) called
#' Distributed Multinomial Regression which utilizes a factorized poisson
#' approximation to the multinomial.  See Roberts, Stewart and Airoldi (2014)
#' for details on the implementation here.  This is dramatically faster than
#' previous versions.  The old default setting which uses a Jeffreys prior is
#' also available.
#' 
#' The argument \code{init.type} allows the user to specify an intialization
#' method. The default 
#' choice, \code{"Spectral"}, provides a deterministic inialization using the
#' spectral algorithm given in Arora et al 2014.  See Roberts, Stewart and
#' Tingley (2016) for details and a comparison of different approaches.
#' Particularly when the number of documents is relatively large we highly
#' recommend the Spectral algorithm which often performs extremely well.  Note
#' that the random seed plays no role in the spectral initialization as it is
#' completely deterministic (unless using the \code{K=0} or random projection
#' settings). When the vocab is larger than 10000 terms we use only the most
#' frequent 10000 terms in creating the initialization.  This may case the 
#' first step of the algorithm to have a very bad value of the objective function
#' but it should quickly stabilize into a good place.  You can tweak the exact 
#' number where this kicks in with the \code{maxV} argument inside control. There
#' appear to be some cases where numerical instability in the Spectral algorithm
#' can cause differences across machines (particularly Windows machines for some reason).
#' It should always give exactly the same answer for a given machine but if you are
#' seeing different answers on different machines, see https://github.com/bstewart/stm/issues/133
#' for a longer explanation.  The other option \code{"LDA"} which uses a few passes
#' of a Gibbs sampler is perfectly reproducible across machines as long as the seed is set.
#' 
#' Specifying an integer greater than 1 for the argument \code{ngroups} causes
#' the corpus to be broken into the specified number of groups.  Global updates
#' are then computed after each group in turn.  This approach, called memoized
#' variational inference in Hughes and Sudderth (2013), can lead to more rapid
#' convergence when the number of documents is large.  Note that the memory
#' requirements scale linearly with the number of groups so this provides a
#' tradeoff between memory efficiency and speed.  The claim of speed here
#' is based on the idea that increasing the number of global updates should
#' help the model find a solution in fewer passes through the document set.
#' However, itt is worth noting that for any particular case the model need 
#' not converge faster and definitely won't converge to the same location. 
#' This functionality should be considered somewhat experimental and we encourage
#'  users to let us know what their experiences are like here in practice.
#' 
#' Models can now be restarted by passing an \code{STM} object to the argument
#' \code{model}.  This is particularly useful if you run a model to the maximum
#' iterations and it terminates without converging.  Note that all the standard
#' arguments still need to be passed to the object (including any formulas, the
#' number of topics, etc.).  Be sure to change the \code{max.em.its} argument
#' or it will simply complete one additional iteration and stop.
#' 
#' You can pass a custom initialization of the beta model parameters to \code{stm}.
#'   
#' 
#' The \code{control} argument is a list with named components which can be
#' used to specify numerous additional computational details.  Valid components
#' include: 
#' \describe{ 
#' \item{\code{tau.maxit}}{Controls the maximum number of
#' iterations when estimating the prior for content covariates.  When the mode
#' is \code{Jeffreys}, estimation proceeds by iterating between the kappa
#' vector corresponding to a particular topic and the associated variance tau
#' before moving on to the next parameter vector. this controls the maximum
#' number of iterations. It defaults to \code{NULL} effectively enforcing
#' convergence.  When the mode is \code{L1} this sets the maximum number of
#' passes in the coordinate descent algorithm and defaults to 1e8.}
#' \item{\code{tau.tol}}{Sets the convergence tolerance in the optimization
#' for content covariates.  When the mode is \code{Jeffreys} this sets the
#' convergence tolerance in the iteration between the kappa vector and
#' variances tau and defaults to 1e-5.  With \code{L1} it defaults to 1e-6.}
#' \item{\code{kappa.mstepmaxit}}{When the mode for content covariate
#' estimation is \code{Jeffreys} this controls the maximum number of passes
#' through the sequence of kappa vectors.  It defaults to 3.  It has no role
#' under \code{L1}- see \code{tau.maxit} option instead.}
#' \item{\code{kappa.msteptol}}{When the mode for content covariate estimation
#' is \code{Jeffreys} this controls the tolerance for convergence (measured by
#' the L1 norm) for the entire M-step.  It is set to .01 by default.  This has
#' no role under mode \code{L1}- see \code{tau.tol} option instead.}
#' \item{\code{fixedintercept}}{a logical indicating whether in content
#' covariate models the intercept should be fixed to the background
#' distribution.  TRUE by default. This only applies when kappa.prior is set to
#' L1.  If FALSE the intercept is estimated from the data without penalty.  In
#' practice estimated intercepts often push term probabilities to zero,
#' resulting in topics that look more like those in a Dirichlet model- that is,
#' most terms have approximately zero probability with some terms with high
#' probability.} 
#' \item{\code{kappa.enet}}{When using the L1 mode for content
#' covariates this controls the elastic net mixing parameter.  See the argument
#' \code{alpha} in \code{glmnet}.  Value must be between 1 and 0 where 1 is the
#' lasso penalty (the default) and 0 is the ridge penalty.  The closer the
#' parameter is to zero the less sparse the solution will tend to be.}
#' \item{\code{gamma.enet}}{Controls the elastic net mixing parameter for the
#' prevalence covariates.  See above for a description.}
#' \item{\code{gamma.ic.k}}{For L1 mode prevalence covariates this controls the 
#' selection of the regularization parameter.  We use a generic information criterion
#'  which penalizes complexity by the parameter \code{ic.k}.  
#'  When set to 2 (as by default) this results in AIC.  When set to log(n) 
#'  (where n is the total number of documents in the corpus) this is equivalent to BIC.
#'    Larger numbers will express a preference for sparser (simpler) models.}
#' \item{\code{gamma.maxits}}{An integer indicating the maximum number of iterations
#' that the prevalence regression variational algorithm can run before erroring out.
#' Defaults to 1000.}
#' \item{\code{nlambda}}{Controls the length of the regularization path when
#' using L1 mode for content covariates.  Defaults to 500.  Note that glmnet
#' relies heavily on warm starts and so a high number will often
#' (counter-intuitively) be less costly than a low number.  We have chosen a
#' higher default here than the default in the glmnet package and we don't
#' recommend changing it.} 
#' \item{\code{lambda.min.ratio}}{For L1 mode content
#' covariates this controls the explored path of regularization values.  This
#' defaults to .0001.  Setting higher numbers will result in more sparse
#' solutions.  This is here primarily for dealing with convergence issues, if
#' you want to favor selection of sparser solutions see the next argument.}
#' \item{\code{ic.k}}{For L1 mode content covariates this controls the
#' selection of the regularization parameter.  We use a generic information
#' criterion which penalizes complexity by the parameter \code{ic.k}.  When set
#' to 2 (as by default) this results in AIC.  When set to log(n) (where n is
#' the total number of words in the corpus) this is equivalent to BIC.  Larger
#' numbers will express a preference for sparser (simpler) models.}
#' \item{\code{nits}}{Sets the number of iterations for collapsed gibbs
#' sampling in LDA initializations.  Defaults to 50} 
#' \item{\code{burnin}}{Sets
#' the burnin for collapsed gibbs sampling in LDA intializations. Defaults to
#' 25} 
#' \item{\code{alpha}}{Sets the prevalence hyperparameter in collapsed
#' gibbs sampling in LDA initializations.  Defaults to 50/K}
#' \item{\code{eta}}{Sets the topic-word hyperparameter in collapsed gibbs
#' sampling in LDa intializations.  Defaults to .01} 
#' \item{\code{contrast}}{A
#' logical indicating whether a standard contrast coding should be used for
#' content covariates.  Typically this should remain at the default of FALSE.}
#' \item{\code{rp.s}}{Parameter between 0 and 1 controlling the sparsity of
#' random projections for the spectral initailization.  Defaults to .05}
#' \item{\code{rp.p}}{Dimensionality of the random projections for the
#' spectral initialization.  Defaults to 3000.}
#' \item{\code{rp.d.group.size}}{Controls the size of blocks considered at a
#' time when computing the random projections for the spectral initialization.
#' Defaults to 2000.} 
#' \item{\code{SpectralRP}}{A logical which when
#' \code{TRUE} turns on the experimental random projections spectral
#' initialization.} 
#' \item{\code{maxV}}{For spectral initializations this will set the maximum
#' number of words to be used in the initialization.  It uses the most frequent words
#' first and then they are reintroduced following initialization.  This allows spectral
#' to be used with a large V.}
#' \item{\code{recoverEG}}{Set to code{TRUE} by default.  If set to \code{FALSE}
#'  will solve the recovery problem in the Spectral algorithm using a downhill simplex
#'  method.  See https://github.com/bstewart/stm/issues/133 for more discussion.}
#' \item{\code{allow.neg.change}}{A logical indicating whether the algorithm is allowed
#' to declare convergence when the change in the bound has become negative. 
#' Defaults to \code{TRUE}.  Set to \code{FALSE} to keep the algorithm from converging
#'  when the bound change is negative.  NB: because this is 
#' only an approximation to the lower-bound the change can be negative at times.  Right
#' now this triggers convergence but the final approximate bound might go higher if you
#' are willing to wait it out. The logic of the default setting is that a negative change
#' in the bound usually means it is barely moving at all.}
#' \item{\code{custom.beta}}{If \code{init.type="Custom"} you can pass your own initialization
#' of the topic-word distributions beta to use as an initialization.  Please note that this takes
#' some care to be sure that it is provided in exactly the right format.  The number of topics and
#' vocab must match exactly.  The vocab must be in the same order.  The values must not be pathological
#' (for instance setting the probability of a single word to be 0 under all topics). The beta should be
#' formatted in the same way as the piece of a returned stm model object \code{stmobj$beta$logbeta}.
#' It should be a list of length the number of levels of the content covariate.  Each element of the list
#' is a K by V matrix containing the logged word probability conditional on the topic.  If you use this
#' option we recommend that you use \code{max.em.its=0} with the model initialization set to random, inspect
#' the returned form of \code{stmobj$beta$logbeta} and ensure that it matches your format.}
#' \item{\code{tSNE_init.dims}}{The K=0 spectral setting uses tSNE to create a low-dimensional
#' projection of the vocab co-occurence matrix.  tSNE starts with a PCA projection as an initialization.
#' We actually do the projection outside the tSNE code so we can use a randomized projection approach.
#' We use the 50 dimensional default of the \pkg{rtsne} package.  That can be changed here.}
#' \item{\code{tSNE_perplexity}}{The \code{Rtsne} function in the \pkg{rtsne} package uses a perplexity
#' parameter.  This defaults to 30 and can throw an error when too high.  \code{stm} will automatically lower
#' the parameter for you until it works, but it can also be directly set here.}
#' }
#' 
#' 
#' @param documents The document term matrix to be modeled. These can be supplied
#' in the native \pkg{stm} format, a sparse term count matrix with one row
#' per document and one column per term, or a
#' \pkg{quanteda} \link[quanteda]{dfm} (document-feature matrix) object.
#' When using the sparse matrix or quanteda format this will include the
#' vocabulary and, for quanteda, optionally the metadata. If using the native list format,
#' the object must be a list of with each element corresponding to a document. Each document is represented
#' as an integer matrix with two rows, and columns equal to the number of unique
#' vocabulary words in the document.  The first row contains the 1-indexed
#' vocabulary entry and the second row contains the number of times that term
#' appears. This is similar to the format in the \code{\link[lda]{lda}} package 
#' except that (following R convention) the vocabulary is indexed from one. Corpora
#' can be imported using the reader function and manipulated using the
#' \code{\link{prepDocuments}}.  Raw texts can be ingested using
#' \code{\link{textProcessor}}. Note that when using \pkg{quanteda} \link[quanteda]{dfm}
#' directly there may be higher memory use (because the texts and metadata are stored
#' twice). You can convert from \pkg{quanteda}'s format directly to our native format
#' using the \pkg{quanteda} function \link[quanteda]{convert}.
#' @param vocab Character vector specifying the words in the corpus in the
#' order of the vocab indices in documents. Each term in the vocabulary index
#' must appear at least once in the documents.  See \code{\link{prepDocuments}}
#' for dropping unused items in the vocabulary.  If \code{documents} is a
#' sparse matrix or \pkg{quanteda} \link[quanteda]{dfm} object, then \code{vocab} should not
#'  (and must not) be supplied.  It is contained already inside the column
#'  names of the matrix.
#' @param K Typically a positive integer (of size 2 or greater) representing
#' the desired number of topics. If \code{init.type="Spectral"} you can also
#' set \code{K=0} to use the algorithm of Lee and Mimno (2014) to set the
#' number of topics (although unlike the standard spectral initialization this
#' is not deterministic).  Additional detail on choosing the number of topics
#' below.
#' @param prevalence A formula object with no response variable or a matrix
#' containing topic prevalence covariates.  Use \code{\link{s}},
#' \code{\link[splines]{ns}} or \code{\link[splines]{bs}} to specify smooth
#' terms. See details for more information.
#' @param content A formula containing a single variable, a factor variable or
#' something which can be coerced to a factor indicating the category of the
#' content variable for each document.
#' @param data an optional data frame containing the prevalence and/or content
#' covariates.  If unspecified the variables are taken from the active
#' environment.
#' @param init.type The method of initialization, by default the spectral initialization.  
#' Must be either Latent
#' Dirichlet Allocation ("LDA"), "Random", "Spectral" or "Custom".  See details for more
#' info. If you want to replicate a previous result, see the argument
#' \code{seed}.  For "Custom" see the format described below under the \code{custom.beta}
#' option of the \code{control} parameters.
#' @param seed Seed for the random number generator. \code{stm} saves the seed
#' it uses on every run so that any result can be exactly reproduced.  When
#' attempting to reproduce a result with that seed, it should be specified
#' here.
#' @param max.em.its The maximum number of EM iterations.  If convergence has
#' not been met at this point, a message will be printed.  If you set this to 
#' 0 it will return the initialization.
#' @param emtol Convergence tolerance.  EM stops when the relative change in
#' the approximate bound drops below this level.  Defaults to .00001.  You 
#' can set it to 0 to have the algorithm run \code{max.em.its} number of steps.
#' See advanced options under \code{control} for more options.
#' @param verbose A logical flag indicating whether information should be
#' printed to the screen.  During the E-step (iteration over documents) a dot
#' will print each time 1\% of the documents are completed.  At the end of each
#' iteration the approximate bound will also be printed.
#' @param reportevery An integer determining the intervals at which labels are
#' printed to the screen during fitting.  Defaults to every 5 iterations.
#' @param LDAbeta a logical that defaults to \code{TRUE} when there are no
#' content covariates.  When set to \code{FALSE} the model performs SAGE style
#' topic updates (sparse deviations from a baseline).
#' @param interactions a logical that defaults to \code{TRUE}.  This
#' automatically includes interactions between content covariates and the
#' latent topics.  Setting it to \code{FALSE} reduces to a model with no
#' interactive effects.
#' @param ngroups Number of groups for memoized inference.  See details below.
#' @param model A prefit model object.  By passing an \code{stm} object to this
#' argument you can restart an existing model.  See details for more info.
#' @param gamma.prior sets the prior estimation method for the prevalence
#' covariate model.  The default \code{Pooled} options uses Normal prior
#' distributions with a topic-level pooled variance which is given a moderately
#' regularizing half-cauchy(1,1) prior.  The alternative \code{L1} uses
#' \code{glmnet} to estimate a grouped penalty between L1-L2.  If your code is running
#' slowly immediately after "Completed E-Step" appears, you may want to switch to the 
#' \code{L1} option. See details below.  
#' @param sigma.prior a scalar between 0 and 1 which defaults to 0.  This sets
#' the strength of regularization towards a diagonalized covariance matrix.
#' Setting the value above 0 can be useful if topics are becoming too highly
#' correlated.
#' @param kappa.prior sets the prior estimation for the content covariate
#' coefficients.  The default option is the \code{L1} prior.  The second option
#' is \code{Jeffreys} which is markedly less computationally efficient but is
#' included for backwards compatability. See details for more information on
#' computation.
#' @param control a list of additional advanced parameters. See details.
#' 
#' @return An object of class STM 
#' 
#' \item{mu}{The corpus mean of topic prevalence and coefficients} 
#' \item{sigma}{Covariance matrix} 
#' \item{beta}{List containing the log of the word probabilities for each topic.}
#' \item{settings}{The settings file. The Seed object will always contain the
#' seed which can be fed as an argument to recover the model.} 
#' \item{vocab}{The vocabulary vector used.} 
#' \item{convergence}{list of convergence elements including the value of the approximate bound on the marginal
#' likelihood at each step.} 
#' \item{theta}{Number of Documents by Number of Topics matrix of topic proportions.} 
#' \item{eta}{Matrix of means for the variational distribution of the multivariate normal latent variables used to
#' calculate theta.} 
#' \item{invsigma}{The inverse of the sigma matrix.}
#' \item{time}{The time elapsed in seconds} 
#' \item{version}{The version number
#' of the package with which the model was estimated.}
#' 
#' @seealso \code{\link{prepDocuments}} \code{\link{labelTopics}}
#' \code{\link{estimateEffect}}
#' @references 
#' Roberts, M., Stewart, B., Tingley, D., and Airoldi, E. (2013)
#' "The structural topic model and applied social science." In Advances in
#' Neural Information Processing Systems Workshop on Topic Models: Computation,
#' Application, and Evaluation. http://goo.gl/uHkXAQ
#' 
#' Roberts M., Stewart, B. and Airoldi, E. (2016) "A model of text for
#' experimentation in the social sciences" Journal of the American Statistical
#' Association.
#' 
#' Roberts, M., Stewart, B., Tingley, D., Lucas, C., Leder-Luis, J., Gadarian,
#' S., Albertson, B., et al. (2014). Structural topic models for open ended
#' survey responses. American Journal of Political Science, 58(4), 1064-1082.
#' http://goo.gl/0x0tHJ
#' 
#' Roberts, M., Stewart, B., & Tingley, D. (2016). "Navigating the Local
#' Modes of Big Data: The Case of Topic Models. In Data Analytics in Social
#' Science, Government, and Industry." New York: Cambridge University Press.
#' @examples
#' 
#' \dontrun{
#' 
#' #An example using the Gadarian data.  From Raw text to fitted model using 
#' #textProcessor() which leverages the tm Package
#' temp<-textProcessor(documents=gadarian$open.ended.response,metadata=gadarian)
#' out <- prepDocuments(temp$documents, temp$vocab, temp$meta)
#' set.seed(02138)
#' mod.out <- stm(out$documents, out$vocab, 3, 
#'                prevalence=~treatment + s(pid_rep), data=out$meta)
#' 
#' #The same example using quanteda instead of tm via textProcessor()
#' #Note this example works with quanteda version 0.9.9-31 and later
#' require(quanteda)
#' gadarian_corpus <- corpus(gadarian, text_field = "open.ended.response")
#' gadarian_dfm <- dfm(gadarian_corpus, 
#'                      remove = stopwords("english"),
#'                      stem = TRUE)
#' stm_from_dfm <- stm(gadarian_dfm, K = 3, prevalence = ~treatment + s(pid_rep),
#'                     data = docvars(gadarian_corpus))
#'                      
#' #An example of restarting a model
#' mod.out <- stm(out$documents, out$vocab, 3, prevalence=~treatment + s(pid_rep), 
#'                data=out$meta, max.em.its=5)
#' mod.out2 <- stm(out$documents, out$vocab, 3, prevalence=~treatment + s(pid_rep), 
#'                 data=out$meta, model=mod.out, max.em.its=10)
#' }
#' @export
stm <- function(documents, vocab, K,
                prevalence=NULL, content=NULL, data=NULL,
                init.type=c("Spectral", "LDA", "Random", "Custom"), seed=NULL,
                max.em.its=500, emtol=1e-5,
                verbose=TRUE, reportevery=5,
                LDAbeta=TRUE, interactions=TRUE,
                ngroups=1, model=NULL,
                gamma.prior=c("Pooled", "L1"), sigma.prior=0,
                kappa.prior=c("L1", "Jeffreys"), control=list())  {
  
  #Match Arguments and save the call
  init.type <- match.arg(init.type)
  Call <- match.call()

  # Convert the corpus to the internal STM format
  args <- asSTMCorpus(documents, vocab, data)
  documents <- args$documents
  vocab <- args$vocab
  data <- args$data
  
  #Documents
  if(missing(documents)) stop("Must include documents")
  if(!is.list(documents)) stop("documents must be a list, see documentation.")
  if(!all(unlist(lapply(documents, is.matrix)))) stop("Each list element in documents must be a matrix. See documentation.")
  if(any(unlist(lapply(documents, function(x) anyDuplicated(x[1,]))))) {
    stop("Duplicate term indices within a document.  See documentation for proper format.")
  }
  N <- length(documents)
  
  #Extract and Check the Word indices
  wcountvec <- unlist(lapply(documents, function(x) rep(x[1,], times=x[2,])),use.names=FALSE)
  #to make this backward compatible we reformulate to old structure.
  wcounts <- list(Group.1=sort(unique(wcountvec)))
  V <- length(wcounts$Group.1)  
  if(!posint(wcounts$Group.1)) {
    stop("Word indices are not positive integers")
  } 
  if(!isTRUE(all.equal(wcounts$Group.1,1:V))) {
    stop("Word indices must be sequential integers starting with 1.")
  } 
  #note we only do the tabulation after making sure it will actually work.
  wcounts$x <- tabulate(wcountvec)
  rm(wcountvec)
  
  #Check the Vocab vector against the observed word indices
  if(length(vocab)!=V) stop("Vocab length does not match observed word indices")
  
  #Check the Number of Topics
  if(missing(K)) stop("K, the number of topics, is required.")
  if(K!=0) {
    #this is the old set of checks
    if(!(posint(K) && length(K)==1 && K>1)) stop("K must be a positive integer greater than 1.")
    if(K==2) warning("K=2 is equivalent to a unidimensional scaling model which you may prefer.")
  } else {
    #this is the special set of checks for Lee and Mimno
    if(init.type!="Spectral") stop("Topic selection method can only be used with init.type='Spectral'")
  }
  #Iterations, Verbose etc.
  if(!(length(max.em.its)==1 & nonnegint(max.em.its))) stop("Max EM iterations must be a single non-negative integer")
  if(!is.logical(verbose)) stop("verbose must be a logical.")
  
  ##
  # A Function for processing prevalence-covariate design matrices
  makeTopMatrix <- function(x, data=NULL) {
    #is it a formula?
    if(inherits(x,"formula")) {
      termobj <- terms(x, data=data)
      if(attr(termobj, "response")==1) stop("Response variables should not be included in prevalence formula.")
      xmat <- try(Matrix::sparse.model.matrix(termobj,data=data),silent=TRUE)
      if(class(xmat)=="try-error") {
        xmat <- try(stats::model.matrix(termobj, data=data), silent=TRUE)
        if(class(xmat)=="try-error") {
                 stop("Error creating model matrix.
                 This could be caused by many things including
                 explicit calls to a namespace within the formula.
                 Try a simpler formula.")
        }
        xmat <- Matrix::Matrix(xmat)
      }
      propSparse <- 1 - Matrix::nnzero(xmat)/length(xmat) 
      #if its less than 50% sparse or there are fewer than 50 columns, just convert to a standard matrix
      if(propSparse < .5 | ncol(xmat) < 50) {
        xmat <- as.matrix(xmat)
      }
      return(xmat)
    }
    if(is.matrix(x)) {
      #Does it have an intercept in first column?
      if(isTRUE(all.equal(x[,1],rep(1,nrow(x))))) return(Matrix::Matrix(x)) 
      else return(cbind(1,Matrix::Matrix(x)))
    }
  }
  
  ###
  #Now we parse both sets of covariates
  ###
  if(!is.null(prevalence)) {
    if(!is.matrix(prevalence) & !inherits(prevalence, "formula")) stop("Prevalence Covariates must be specified as a model matrix or as a formula")
    xmat <- makeTopMatrix(prevalence,data)
    if(is.na(nnzero(xmat))) stop("Missing values in prevalence covariates.")
  } else {
    xmat <- NULL
  }
  
  if(!is.null(content)) {
    if(inherits(content, "formula")) {
      termobj <- terms(content, data=data)
      if(attr(termobj, "response")==1) stop("Response variables should not be included in content formula.")
      if(nrow(attr(termobj, "factors"))!=1) stop("Currently content can only contain one variable.")
      if(is.null(data)) {
        yvar <- eval(attr(termobj, "variables"))[[1]]
      } else {
        char <- rownames(attr(termobj, "factors"))[1]
        yvar <- data[[char]]
      }
      yvar <- as.factor(yvar)
    } else {
      yvar <- as.factor(content)
    }
    if(any(is.na(yvar))) stop("Your content covariate contains missing values.  All values of the content covariate must be observed.")
    yvarlevels <- levels(yvar)
    betaindex <- as.numeric(yvar)
  } else{
    yvarlevels <- NULL
    betaindex <- rep(1, length(documents))
  }
  A <- length(unique(betaindex)) #define the number of aspects
  
  #Checks for Dimension agreement
  ny <- length(betaindex)
  nx <- ifelse(is.null(xmat), N, nrow(xmat))
  if(N!=nx | N!=ny) stop(paste("number of observations in content covariate (",ny,
                               ") prevalence covariate (",
                               nx,") and documents (",N,") are not all equal.",sep=""))
  
  #Some additional sanity checks
  if(!is.logical(LDAbeta)) stop("LDAbeta must be logical")
  if(!is.logical(interactions)) stop("Interactions variable must be logical")
  if(sigma.prior < 0 | sigma.prior > 1) stop("sigma.prior must be between 0 and 1")
  if(!is.null(model)) {
    if(max.em.its <= model$convergence$its) stop("when restarting a model, max.em.its represents the total iterations of the model 
                                                 and thus must be greater than the length of the original run")
  }
  ###
  # Now Construct the Settings File
  ###
  settings <- list(dim=list(K=K, A=A, 
                            V=V, N=N, wcounts=wcounts),
                   verbose=verbose,
                   topicreportevery=reportevery,
                   convergence=list(max.em.its=max.em.its, em.converge.thresh=emtol, 
                                    allow.neg.change=TRUE),
                   covariates=list(X=xmat, betaindex=betaindex, yvarlevels=yvarlevels, formula=prevalence),
                   gamma=list(mode=match.arg(gamma.prior), prior=NULL, enet=1, ic.k=2,
                              maxits=1000),
                   sigma=list(prior=sigma.prior),
                   kappa=list(LDAbeta=LDAbeta, interactions=interactions, 
                              fixedintercept=TRUE, mstep=list(tol=.001, maxit=3),
                              contrast=FALSE),
                   tau=list(mode=match.arg(kappa.prior), tol=1e-5,
                            enet=1,nlambda=250, lambda.min.ratio=.001, ic.k=2,
                            maxit=1e4),
                   init=list(mode=init.type, nits=50, burnin=25, alpha=(50/K), eta=.01,
                             s=.05, p=3000, d.group.size=2000, recoverEG=TRUE,
                             tSNE_init.dims=50, tSNE_perplexity=30), 
                   seed=seed,
                   ngroups=ngroups)
  if(init.type=="Spectral" & V > 10000) {
    settings$init$maxV <- 10000
  }
  
  if(settings$gamma$mode=="L1") {
    #if(!require(glmnet) | !require(Matrix)) stop("To use L1 penalization please install glmnet and Matrix")
    if(ncol(xmat)<=2) stop("Cannot use L1 penalization in prevalence model with 2 or fewer covariates.")
  }

    
  ###
  # Fill in some implied arguments.
  ###
  
  #Is there a covariate on top?
  if(is.null(prevalence)) {
    settings$gamma$mode <- "CTM" #without covariates has to be estimating the mean.
  } 
  
  #Is there a covariate on the bottom?
  if(is.null(content)) {
    settings$kappa$interactions <- FALSE #can't have interactions without a covariate.
  } else {
    settings$kappa$LDAbeta <- FALSE #can't do LDA topics with a covariate 
  }
  
  ###
  # process arguments in control
  ###
  
  #Full List of legal extra arguments
  legalargs <-  c("tau.maxit", "tau.tol", 
                  "fixedintercept","kappa.mstepmaxit", "kappa.msteptol", 
                  "kappa.enet", "nlambda", "lambda.min.ratio", "ic.k", "gamma.enet",
                  "gamma.ic.k",
                  "nits", "burnin", "alpha", "eta", "contrast",
                  "rp.s", "rp.p", "rp.d.group.size", "SpectralRP",
                  "recoverEG", "maxV", "gamma.maxits", "allow.neg.change",
                  "custom.beta", "tSNE_init.dims", "tSNE_perplexity")
  if (length(control)) {
    indx <- pmatch(names(control), legalargs, nomatch=0L)
    if (any(indx==0L))
      stop(gettextf("Argument %s not matched", names(control)[indx==0L]),
           domain = NA)
    fullnames <- legalargs[indx]
    for(i in fullnames) {
      if(i=="tau.maxit") settings$tau$maxit <- control[[i]]
      if(i=="tau.tol") settings$tau$tol <- control[[i]]
      if(i=="fixedintercept")settings$kappa$fixedintercept <- control[[i]]
      if(i=="kappa.enet") settings$tau$enet <- control[[i]]
      if(i=="kappa.mstepmaxit") settings$kappa$mstep$maxit <- control[[i]] 
      if(i=="kappa.msteptol") settings$kappa$mstep$tol <- control[[i]] 
      if(i=="nlambda") settings$tau$nlambda <- control[[i]]
      if(i=="lambda.min.ratio") settings$tau$lambda.min.ratio <- control[[i]]
      if(i=="ic.k") settings$tau$ic.k <- control[[i]]
      if(i=="gamma.enet") settings$gamma$enet <- control[[i]]
      if(i=="gamma.ic.k") settings$gamma$ic.k <- control[[i]]
      if(i=="nits") settings$init$nits <- control[[i]]
      if(i=="burnin") settings$init$burnin <- control[[i]]
      if(i=="alpha") settings$init$alpha <- control[[i]]
      if(i=="eta") settings$init$eta <- control[[i]]
      if(i=="contrast") settings$kappa$contrast <- control[[i]]
      if(i=="rp.s")  settings$init$s <- control[[i]]
      if(i=="rp.p")  settings$init$p <- control[[i]]
      if(i=="rp.d.group.size")  settings$init$d.group.size <- control[[i]]
      if(i=="SpectralRP" && control[[i]]) settings$init$mode <- "SpectralRP" #override to allow spectral rp mode
      if(i=="recoverEG" && !control[[i]]) settings$init$recoverEG <- control[[i]]
      if(i=="maxV" && control[[i]]) {
        settings$init$maxV <- control[[i]]
        if(settings$init$maxV > V) stop("maxV cannot be larger than the vocabulary")
      }
      if(i=="tSNE_init.dims" && control[[i]]) settings$init$tSNE_init.dims <- control[[i]]
      if(i=="tSNE_perplexity" && control[[i]]) settings$init$tSNE_perplexity <- control[[i]]
      if(i=="gamma.maxits") settings$gamma$maxits <- control[[i]]
      if(i=="allow.neg.change") settings$convergence$allow.neg.change <- control[[i]]
      if(i=="custom.beta") {
        if(settings$init$mode!="Custom") {
          warning("Custom beta supplied, setting init argument to Custom.")
          settings$init$mode <- "Custom"
        }
        settings$init$custom <- control[[i]]
      }
    }
  }
  
  ###
  # Process the Seed
  ###
  if(is.null(settings$seed)) {
    #if there is no seed, choose one and set it, recording for later
    seed <- floor(runif(1)*1e7) 
    set.seed(seed)
    settings$seed <- seed
  } else {
    #otherwise just use the provided seed.
    set.seed(settings$seed)
  }
  
  settings$call <- Call
  ###
  # Finally run the actual model
  ###
  return(stm.control(documents, vocab, settings,model))
}
