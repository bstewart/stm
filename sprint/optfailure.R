#Load up the data
load("optfailure.RData")
#contains 6 existing model objects containing 2 trajectories
#of the same model.
# and their documents

## Remember optim is minimizing so a lower objective is better!

#define a function that fits a given document 
#using model parameters while allowing control
#of outputs and various inputs
fitme <- function(dn, mod, documents,
                  control=list(maxit=500), 
                  method="BFGS",
                  alt_init=NULL,hess=FALSE, eval_only=FALSE) {
  ### Arguments
  # dn          document number
  # mod         stm model object
  # documents   documents
  # control     control arguments passed to optim
  # method      the optimization argument for optim.
  # alt_init    an alternative initialization vector
  # hess        if true calculate the hessian
  # eval_only   no optimization, just evaluates at the
  #             initialization. if hess=TRUE will also
  #             calculate the gradient and negative hessian.
  ###
  doc <- documents[[dn]]
  words <- doc[1,]
  aspect <- mod$settings$covariates$betaindex[dn]
  if(is.null(alt_init)) {
    init <- mod$eta[dn,]
  } else {
    init <- alt_init
  }
  mu.i <- mod$mu$mu[,dn]
  beta.i <- exp(mod$beta$logbeta[[aspect]][,words,drop=FALSE])
  sigobj <- try(chol.default(mod$sigma), silent=TRUE)
  if(class(sigobj)=="try-error") {
    sigmaentropy <- (.5*determinant(sigma, logarithm=TRUE)$modulus[1])
    siginv <- solve(sigma)
  } else {
    sigmaentropy <- sum(log(diag(sigobj)))
    siginv <- chol2inv(sigobj)
  }
  #if we only want evaluation just return initialization
  if(eval_only) {
    if(hess) {
      grad <- stm:::gradcpp(eta = init, beta = beta.i, doc_ct = doc[2,],
                            mu=mu.i, siginv=siginv)
      hess <- stm:::hpbcpp2(eta = init, beta = beta.i, doc_ct = doc[2,],
                            mu=mu.i, siginv=siginv,sigmaentropy=sigmaentropy)
    } else {
      hess <- NULL
      grad <- NULL
    }
    return(list(param=init, 
                obj=stm:::lhoodcpp(init, beta.i, doc[2,], mu.i, siginv),
                grad=grad,
                hess=hess[[1]]))
  }
  
  #infer the document
  doc.results <- stm:::logisticnormalcpp(eta=init, mu=mu.i, siginv=siginv, beta=beta.i, 
                                         doc=doc, sigmaentropy=sigmaentropy, control=control,
                                         method=method,hpbcpp=hess)
  list(param=doc.results$eta$lambda, 
       obj=stm:::lhoodcpp(doc.results$eta$lambda, beta.i, doc[2,], mu.i, siginv),
       hess=doc.results)
}
#a softmax function
softmax0 <- function(x) stm:::softmax(c(x,0))

#Here we have 2 models, both started from the same initialization
#they both use the Neumaier sum for sigma and beta updating and
#randomize order of document updating.  They diverge.  The two
#models are "a" and "b" with version separated by an underscore
#indicating the iteration where they were saved.  Thus "a_1" is
#model "a" after 1 iteration, and "b_3" is model "b" after 3 
#iteration.  Again the only difference between "a" and "b" is the
#random seed determining document order.

##
#Let's check for divergence in theta
##
#after 1 iteration, it is undetectable at machine precision
all.equal(a_1$theta,b_1$theta,tol=.Machine$double.eps)
#after 2nd iteration, it is small 4.25e-08
all.equal(a_2$theta,b_2$theta)
#after 3rd iteration, there is a noticeable difference!
all.equal(a_3$theta,b_3$theta)

#the difference in third iteration is almost .2 on a single
#document loading!  Yikes!
max(abs(a_3$theta-b_3$theta))
#this is happening on document 11,647
which(abs(a_3$theta-b_3$theta)>.1,arr.ind=TRUE)
#let's try to understand how this happens...

#Third iteration deviations come from 2nd iteration global
#parameters. I made a function to refit the local parameters
#of a document based on a given model.
a2doc <- fitme(11647, a_2, documents)
b2doc <- fitme(11647, b_2, documents)
#differences are concentrated in topics 1,15,99 and they are big:
softmax0(a2doc$par)[c(1,15,99)]
softmax0(b2doc$par)[c(1,15,99)]
#the difference is in the initialization, not the model parameters
#initializing under b2 solution gives a bad answer for a2, initializing
#under the a2 solution gives a better answer for b2 parameters.
##a2 and b2 parameters with b2 init and b2 solution init
fitme(11647, a_2, documents,alt_init = b2doc$param)$obj #931.911
fitme(11647, b_2, documents)$obj #931.911
##a2 and b2 parameters both have better solutions available:
fitme(11647, a_2, documents)$obj #930.795
fitme(11647, b_2, documents,alt_init = a2doc$param)$obj #930.795



#there are 4 possible inputs that could cause the differences here
# - sigma
# - beta
# - eta (the initialization)
# - mu

#Let's see if they are different:
all.equal(a_2$sigma,b_2$sigma) #similar at e-7
all.equal(a_2$beta, b_2$beta, tol=.Machine$double.eps) #similar at e-9
all.equal(a_2$eta[11647,],b_2$eta[11647,], tol=.Machine$double.eps) #similar at e-16
all.equal(a_2$mu$mu[,11647],b_2$mu$mu[,11647]) #similar at e-7
#they are very similar and I had to drop to machine tolerance to 
#see any difference in some of them.

#How do we get such different results, then?  Let's try to see what
#inputs have to change to get different answers
ablate <- function(vec, orig_mod, new_mod,documents,
                   docnum=11647,topics=c(1,15,99)) {
  #this function takes an original model and depending
  #on its arguments changes pieces of the global parameters
  #by replacing them with those from newmodel.
  #it then outputs results for specific topics and the objective function
  #note that eta is the initialization
  check <- orig_mod
  if(vec[1]>0)  check$sigma <- new_mod$sigma
  if(vec[2]>0) check$beta <- new_mod$beta
  if(vec[3]>0) check$eta <- new_mod$eta
  if(vec[4]>0) check$mu <- new_mod$mu
  m <- fitme(docnum,check,documents)
  c(softmax0(m$param)[topics],m$obj)
}
#this next line creates every possible combination of inputs
val <- expand.grid(c(0,1),c(0,1),0:1,0:1)
colnames(val) <- c("sigma","beta","init","mu")
ablation <- t(apply(val, 1, ablate, orig_mod=a_2,new_mod=b_2, documents=documents))
ablation #essentially the same but lines 1 and 5
#looking at this we can see essentially the only thing that
#doesn't get the answer we get with b_2 is if we change the
#initialization.
# In other words: any mild change in the parameters gets you a worse
# solution.  The only thing you can change is the initialization which
# if you don't change any parameters, will still give you the better solution
ablation <- t(apply(val, 1, ablate, orig_mod=b_2,new_mod=a_2, documents=documents))
ablation #essentially the same but lines 1 and 5

#it turns out actually, you can make incredibly subtle changes
check <- a_2
#look at these values for mu
check$mu$mu[c(31,95),11647]
#nudge them up by 1e-6
check$mu$mu[c(31,95),11647] <- check$mu$mu[c(31,95),11647] + c(1e-6,1e-6)
#and the result is completely different
softmax0(fitme(11647,check,documents)$param)[c(1,15,99)]
#vs.
softmax0(fitme(11647,a_2,documents)$param)[c(1,15,99)]



####
# Optimization failure
####
#This is an optimization failure.
# we can see that because we can
# return a2's better solution by changing
# the initialization of "check" to the final
# a2 parameters

fitme(11647,a_2,documents)$obj #930.7915
fitme(11647,check,documents)$obj #931.911
fitme(11647,check,documents, 
      alt_init = fitme(11647,a_2,documents)$param)$obj #930.7915
#note that changing to a2's own initialization is 
#NOT sufficient to get the better answer with the parameter
#changes. It has to be the better solution.
fitme(11647,check,documents,
      alt_init=a_2$eta[11647,])$obj #931.911

#define a quick function to spit out the objective and the
#softmax for the three problem topics
look <- function(x,topics=c(1,15,99)) list(obj=x$obj,par=softmax0(x$param)[topics])

##
# BFGS seemingly behaves different from the other optimizers
#alternate optimization
look(fitme(11647,a_2,documents,method="L-BFGS-B")) #931.911
look(fitme(11647,a_2,documents,method="BFGS")) #930.7915 (best!)
look(fitme(11647,a_2,documents,method="CG")) #931.911 
look(fitme(11647,a_2,documents,method="Nelder-Mead")) #934.5557
look(fitme(11647,a_2,documents,method="SANN")) #934.5874
# I have absolutely no explanation for why

#As best I can tell, this isn't an issue of convergence tolerances that are too loose
look(fitme(11647,a_2,documents,method="CG", control=list(reltol=.Machine$double.eps,
                                                         maxit=5000))) 

#an interesting open question here:
#   The optimization failures are triggered by very small differences.
#   Why?
#   I would hypothesize from this that the differences are very, very 
#  small but get magnified by these big optimization problems.
#
#  when calculating the Hessian (see STMCfuns.cpp) if it isn't positive
#  definite, we do diagonal dominance to force it to be invertible. Would
#  checking this identify these problems?
#  can we fix these problems with a better optimization routine in general?
#  if not, could we trigger an alterantive routine that only works when
#  the hessian isn't invertible?



###
# Below here I am just messing around
###

x <- fitme(11647,a_2,documents, eval_only=TRUE,hess=TRUE)
y <- fitme(11647,check,documents,eval_only=TRUE,hess=TRUE)
x$obj; y$obj
x$grad -y$grad #there are barely any appreciable differences here...
## IDEA: test the gradient path and the hessian*gradient path from
# these start points and plot out the objective...


#best known parameters we have
best <- fitme(11647,a_2,documents,method="BFGS")$par
bad <- fitme(11647,check,documents,method="BFGS")$par
#best objective
fitme(11647,a_2,documents,method="BFGS", alt_init=best, eval_only=TRUE)$obj
#if you give it the bad one it doesn't move
fitme(11647,a_2,documents,method="BFGS", alt_init=bad, eval_only=TRUE)$obj
fitme(11647,a_2,documents,method="BFGS", alt_init=bad)$obj

diff <- best - bad
c <- seq(-.2,1.2,by=.01)
obj <- vector(mode="numeric",length=length(c))
for(i in 1:length(c)) {
  obj[i] <- fitme(11647,a_2,documents, alt_init=bad + c[i]*diff, eval_only=TRUE)$obj
}
plot(c,obj)


x <- fitme(11647,a_2,documents,alt_init=bad, eval_only=TRUE,hess=TRUE)
y <- fitme(11647,a_2,documents,alt_init=best, eval_only=TRUE,hess=TRUE)
dir <- -chol2inv(chol(x$hess))%*%x$grad
x_a <- fitme(11647,a_2,documents,alt_init=alt, eval_only=TRUE,hess=TRUE)

newton <- function(init,n) {
  obj <- vector(mode="numeric",length=n)
  opts <- vector(mode="numeric",length=n)
  par <- init
  for(i in 1:n) {
    x <- fitme(11647,a_2,documents,alt_init=par, eval_only=TRUE,hess=TRUE)
    obj[i] <- x$obj 
    if(i==n) return(list(obj=obj,par=par,opts=opts))
    dir <- -chol2inv(chol(x$hess))%*%x$grad
    f <- function(x) stm:::lhoodcpp(par + x*dir, beta.i, doc[2,], mu.i, siginv)
    c <- optimize(f, interval=c(-10,10))
    opts[i] <- c$minimum
    par <- par + c$minimum*dir
  }
}
better <- newton(bad,5000)  
better$opts
better$obj[10000] - better$obj[9999]
better$obj[5000] - better$obj[1]

x_step <- fitme(11647,a_2,documents,alt_init=bad - dir, eval_only=TRUE,hess=TRUE)

exp(a_2$beta$logbeta[[1]])[1,]

stm:::labelTopics(a_2, topics=c(1,99))

diag(chol2inv(chol(y$hess)))[c(1,15,99)]


#can we consider small deviations?
diff <- best - bad
print(fitme(11647,a_2,documents,method="BFGS", alt_init=bad + diff*.1)$obj,10)
##check hessian? look at gradient calculation?