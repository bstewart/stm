##########################
# CONSTANTS
##########################
sigma.round.each = 8
beta.round.each = 8
sigma.round = 8
beta.round = 8

control = list(estep.sigma.round=sigma.round, estep.sigma.round.each=sigma.round.each, estep.beta.round=beta.round, estep.beta.round.each=beta.round.each)
##########################

library(stm)

seed <- sample(1:1000000, 1)
# seed <- 517966
cat('------------------------')
cat(seed)
cat('------------------------')

m1 <- stm(poliblog5k.docs, poliblog5k.voc, K=10, seed=seed, cores=1, control=control)
m2 <- stm(poliblog5k.docs, poliblog5k.voc, K=10, seed=seed, cores=7, control=control)
all.equal(m1, m2)
