library(stm)

deter <- stm(poliblog5k.docs, poliblog5k.voc, K=20, prevalence=~rating, data=poliblog5k.meta, max.em.its = 100)

rand <- stm(poliblog5k.docs, poliblog5k.voc, K=20, prevalence=~rating, data=poliblog5k.meta, seed=13, max.em.its = 100, control=list(randomize=TRUE))

all.equal(deter, rand)

