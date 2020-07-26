#!/usr/bin/env Rscript

# Get command line argument for operating system
args = commandArgs(trailingOnly=TRUE)

file_names = c("../artifact/reg_Linux.Rds", "../artifact/reg_macOS.Rds", "../artifact/reg_Windows.Rds", "../artifact/neum_Linux.Rds", "../artifact/neum_macOS.Rds", "../artifact/neum_Windows.Rds")

for(i in 1:length(file_names)) {
  if(file_names[i] == "../artifact/reg_Linux.Rds") {
    reg_Linux <- readRDS(file_names[i])
  }
  else if(file_names[i] == "../artifact/neum_macOS.Rds") {
    neum_macOS <- readRDS(file_names[i])
  }
}

print(all.equal(reg_Linux, neum_macOS, tolerance=1.0e-13))
