library(stm)

num_topics <- numeric()
avg_time_sum <- numeric()
num_runs <- numeric()
sum_type <- character()

ks <- c(5, 10, 20, 40, 80, 100, 120)
num_run = 3

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
  sum_type <- cbind(sum_type, "reg")
}

reg_sum_df <- data.frame(
  num_topics = as.numeric(num_topics),
  avg_time_sum = as.numeric(avg_time_sum),
  num_runs = as.numeric(num_runs),
  sum_type = as.character(sum_type)
)

#print(reg_sum_df)
save(reg_sum_df,file="/home/gyorgym/sociology/stm/sprint/reg_sum.Rda")
