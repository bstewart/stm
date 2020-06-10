library(stm)

# Get the docs
docs = poliblog5k.docs
voc = poliblog5k.voc
heldout <- make.heldout(docs, voc)

# For dataframe
num_topics <- numeric()
num_iter <- numeric()
var_bound <- numeric()
doc_order <- character()

ks <- 5:100

cat("Summation time trials\n")
for(k in ks) {
  cat(paste("k:", k, "\n"))
  cat("deterministic run\n")
  deter <- stm(heldout$documents, heldout$vocab, K=k, init.type = "Spectral", max.em.its = 100, verbose = FALSE, control=list(method="BFGS", neum_sum_cpp=TRUE)) 
  num_topics <- cbind(num_topics, k)
  num_iter <- cbind(num_iter, deter$convergence$its)
  var_bound <- cbind(var_bound, max(deter$convergence$bound))
  doc_order <- cbind(doc_order, "deter")    
  cat("randomized run\n")
  rand <- stm(heldout$documents, heldout$vocab, K=k, init.type = "Spectral", seed=13, max.em.its = 100, verbose = FALSE, control=list(method="BFGS", neum_sum_cpp=TRUE, rand_docs=TRUE)) 
  num_topics <- cbind(num_topics, k)
  num_iter <- cbind(num_iter, rand$convergence$its)
  var_bound <- cbind(var_bound, max(rand$convergence$bound))
  doc_order <- cbind(doc_order, "rand")    
  print(all.equal(deter$phi, rand$phi))
  print(all.equal(deter$sigma, rand$sigma))
  print(all.equal(deter$theta, rand$theta))
  print(all.equal(deter$nu, rand$nu))
}

neum_sum_df <- data.frame(
  num_topics = as.numeric(num_topics),
  num_iter = as.numeric(num_iter),
  var_bound = as.numeric(var_bound),
  doc_order = as.character(sum_type)
)

save(neum_sum_df,file="neum_sum_df.Rda")