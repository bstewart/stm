library(stm)

num_topics <- numeric()
avg_time_sum <- numeric()
num_runs <- numeric()
num_iter <- numeric()
var_bound <- numeric()
sum_type <- character()

ks <- c(40) #5, 10, 20, 40, 80)
num_run = 1

print("Summation time trials")
for(k in ks) {
  cum_time = 0
  for(j in 1:num_run) {
    cat(paste("k:", k, "j:", j, "\n"))
    res <- stm(poliblog5k.docs, poliblog5k.voc, K=k, prevalence=~rating, data=poliblog5k.meta, max.em.its = 100, verbose = FALSE, control=list(method="BFGS"))
    cum_time = cum_time + res$time
  }
  num_topics <- cbind(num_topics, k)
  avg_time_sum <- cbind(avg_time_sum, cum_time/num_run)
  num_runs <- cbind(num_runs, num_run)
  num_iter <- cbind(num_iter, res$convergence$its)
  var_bound <- cbind(var_bound, max(res$convergence$bound))
  sum_type <- cbind(sum_type, "reg")
}

reg_sum_df <- data.frame(
  num_topics = as.numeric(num_topics),
  avg_time_sum = as.numeric(avg_time_sum),
  num_iter = as.numeric(num_iter),
  var_bound = as.numeric(var_bound),
  num_runs = as.numeric(num_runs),
  sum_type = as.character(sum_type)
)

#print(reg_sum_df)
save(reg_sum_df,file="/home/gomfy/sociology/stm/git/sprint/stm/sprint/reg_sum.Rda")
