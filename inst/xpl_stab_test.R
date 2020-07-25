#!/usr/bin/env Rscript
library(stm)
library(stringr)

# Get command line argument for operating system
args = commandArgs(trailingOnly=TRUE)

file_names = c("reg_oneLinux.Rds", "reg_onemacOS.Rds", "reg_oneWindows.Rds", "neum_oneLinux.Rds", "neum_onemacOS.Rds", "neum_oneWindows.Rds")

for(i in 1:length(file_names)) {
  if(file_names[i] == "reg_oneLinux.Rds") {
    reg_oneLinux <- readRDS(file_names[i])
  }
  else if(file_names[i] == "neum_onemacOS.Rds") {
    neum_onemacOS <- readRDS(file_names[i])
  }
}

print(all.equal(reg_oneLinux, neum_onemacOS, tolerance=1.0e-13))
