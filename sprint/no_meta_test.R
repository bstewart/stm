library(stm)

heldout <- make.heldout(poliblog5k.docs, poliblog5k.voc)

deter <- stm(heldout$documents, heldout$vocab, K = 5, init.type = "Spectral", max.em.its = 100) 
rand <- stm(heldout$documents, heldout$vocab, K = 5, init.type = "Spectral", seed=13, max.em.its = 100, control=list(randomize=TRUE))

all.equal(deter, rand)
