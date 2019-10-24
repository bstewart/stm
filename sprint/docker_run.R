library(stringr)
library(stm)

mod <- stm(poliblog5k.docs, poliblog5k.voc, K=10, prevalence=~rating, data=poliblog5k.meta, max.em.its = 100)

fname = str_c("../output/",format(Sys.time(), "%H-%M-%S"), ".mod")
save(mod, file = fname)

