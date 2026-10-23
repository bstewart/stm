#!/usr/bin/env Rscript

source(file.path("benchmarks", "large-scale-common.R"))
options <- parse_benchmark_args(list(
  results="large-scale-results",
  output="large-scale-summary.csv"
))

files <- list.files(options$results, pattern="[.]rds$", full.names=TRUE)
if(!length(files)) stop("No RDS results found in ", options$results)
rows <- lapply(files, function(path) {
  result <- readRDS(path)
  if(is.null(result$summary)) return(NULL)
  result$summary$result.file <- basename(path)
  result$summary
})
summary <- do.call(rbind, rows)
summary <- summary[order(summary$seed, summary$method),,drop=FALSE]
utils::write.csv(summary, options$output, row.names=FALSE)
print(summary)

