library(stm)

detercpphm <- stm(poliblog5k.docs, poliblog5k.voc, K=5, prevalence=~rating, data=poliblog5k.meta, max.em.its = 100, control=list(method="BFGS", neum_sum_cpp=TRUE))

#deter <- stm(poliblog5k.docs, poliblog5k.voc, K=5, prevalence=~rating, data=poliblog5k.meta, max.em.its = 100, control=list(method="BFGS"))
#randcpphm <- stm(poliblog5k.docs, poliblog5k.voc, K=5, prevalence=~rating, data=poliblog5k.meta, seed=13, max.em.its = 100, control=list(method="BFGS", rand_docs=TRUE, neum_sum_cpp=TRUE))

all.equal(detercpphm, randcpp)

