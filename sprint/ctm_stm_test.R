library(topicmodels)
library(stm)

K <- 6
mag <- 2138
set.seed(2138)
heldout <- make.heldout(poliblog5k.docs, poliblog5k.voc)
slam <- convertCorpus(heldout$documents, heldout$vocab, type="slam")

ctm <- CTM(slam, k = K, control=list(seed=as.integer(2138),verbose=1))
deter <- stm(heldout$documents, heldout$vocab, K = K, init.type = "Spectral", max.em.its = 100)
rand <- stm(heldout$documents, heldout$vocab, K = K, init.type = "Spectral", max.em.its = 100, control=list(randomize=TRUE))

all.equal(deter, rand)
