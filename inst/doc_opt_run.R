#!/usr/bin/env Rscript
library(stm)
library(stringr)

# Get command line argument for operating system
args = commandArgs(trailingOnly=TRUE)

mod <- readRDS("./inst/mod.Rds")
method <- "ucminf"
tol <- 1e-6
file_name = str_c(method, "_1e-6_", args[1], ".Rds")
#file_name = str_c("./inst/docopt/", method, "_1e-6_", "Windows", ".Rds")

doc_id <- numeric()
siginv_norm2 <- numeric()
obj_value <- numeric()
sol_norm1 <- numeric()
sol_norm2 <- numeric()

# Optimize every document and save solution stats
for(i in 1:length(poliblog5k.docs)) {
    docopt <- stm:::optimize_doc(dn = i, 
                                 mod = mod, 
                                 documents = poliblog5k.docs, 
                                 method = method,
                                 tol = tol)
    doc_id <- cbind(doc_id, i)
    siginv_norm2 <- cbind(siginv_norm2, norm(as.matrix(docopt$siginv), type = "f"))
    obj_value <- cbind(obj_value, docopt$optim.out$value)
    sol_norm1 <- cbind(sol_norm1, max(abs(docopt$optim.out$par)))
    sol_norm2 <- cbind(sol_norm2, norm(as.matrix(docopt$optim.out$par), type = "f"))
    # cat(paste("value:", docopt$optim.out$value, "(tol:", tol, ")\n"))
    # cat("solution: \n")
    # print(docopt$optim.out$par)
    # cat(paste("norm 1: ", max(abs(docopt$optim.out$par)), "\n"))
    # cat(paste("norm 2: ", norm(as.matrix(docopt$optim.out$par), type = "f"), "\n"))
}

doc_opt_df <- data.frame(
  doc_id = as.numeric(doc_id),
  siginv_norm2 = as.numeric(siginv_norm2),
  obj_value = as.numeric(obj_value),
  sol_norm1 = as.numeric(sol_norm1),
  sol_norm2 = as.numeric(sol_norm2)
)

#print(doc_opt_df)
saveRDS(doc_opt_df, file=file_name)


