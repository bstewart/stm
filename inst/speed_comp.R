library(stm)

# Get the docs
docs = poliblog5k.docs
voc = poliblog5k.voc
heldout <- make.heldout(docs, voc)

# For dataframe
num_topics <- numeric()
num_iter <- numeric()
run_time <- numeric()
var_bound <- numeric()
sum_type <- character()

ks <- 5:100

cat("Summation speed comparisons\n")
for(k in ks) {
  cat(paste("k:", k, "\n"))
  cat("Determine max.em.its\n")
  neum <- stm(heldout$documents, heldout$vocab, K=k, init.type = "Spectral", max.em.its = 100, verbose = FALSE, control=list(method="BFGS", neum_sum_cpp=TRUE)) 
  neum_iter <- neum$convergence$its
  reg <- stm(heldout$documents, heldout$vocab, K=k, init.type = "Spectral", max.em.its = 100, verbose = FALSE, control=list(method="BFGS")) 
  reg_iter <- reg$convergence$its
  max_iter <- min(neum_iter, reg_iter)
  
  cat("Neumaier sum run\n")
  neum <- stm(heldout$documents, heldout$vocab, K=k, init.type = "Spectral", max.em.its = max_iter, verbose = FALSE, control=list(method="BFGS", neum_sum_cpp=TRUE)) 
  num_topics <- cbind(num_topics, k)
  num_iter <- cbind(num_iter, neum$convergence$its)
  run_time <- cbind(run_time, neum$time)
  var_bound <- cbind(var_bound, max(neum$convergence$bound))
  sum_type <- cbind(sum_type, "neum")    
  
  cat("Regular sum run\n")
  reg <- stm(heldout$documents, heldout$vocab, K=k, init.type = "Spectral", max.em.its = max_iter, verbose = FALSE, control=list(method="BFGS")) 
  num_topics <- cbind(num_topics, k)
  num_iter <- cbind(num_iter, reg$convergence$its)
  run_time <- cbind(run_time, reg$time)
  var_bound <- cbind(var_bound, max(reg$convergence$bound))
  sum_type <- cbind(sum_type, "reg")    
}

speed_comp_df <- data.frame(
  num_topics = as.numeric(num_topics),
  num_iter = as.numeric(num_iter),
  run_time = as.numeric(run_time),
  var_bound = as.numeric(var_bound),
  sum_type = as.character(sum_type)
)

save(speed_comp_df,file="speed_comp_df.Rda")