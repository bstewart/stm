#Simulates betas for each topic
#Takes parameters value from estimateEffect
#nsims is the number of simulations needed

simBetas <- function(parameters, nsims=100){
  simbetas <- list()
  for(i in 1:length(parameters)) {
      simbetas[[i]] <- do.call(rbind, lapply(parameters[[i]],
                                             function(x) rmvnorm(n=nsims, mu=x$est, Sigma=x$vcov)))
    }
  return(simbetas)
}
