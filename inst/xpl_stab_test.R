#!/usr/bin/env Rscript

# Get command line argument for operating system
args = commandArgs(trailingOnly=TRUE)

file_names = c("./artifact/reg_Linux.Rds", "./artifact/reg_macOS.Rds", "./artifact/reg_Windows.Rds", "./artifact/neum_Linux.Rds", "./artifact/neum_macOS.Rds", "./artifact/neum_Windows.Rds")

for(i in 1:length(file_names)) {
  if(file_names[i] == "./artifact/reg_Linux.Rds") {
    reg_Linux <- readRDS(file_names[i])
  }
  else if(file_names[i] == "./artifact/reg_macOS.Rds") {
    reg_macOS <- readRDS(file_names[i])
  }
  else if(file_names[i] == "./artifact/reg_Windows.Rds") {
    reg_Windows <- readRDS(file_names[i])
  }
  else if(file_names[i] == "./artifact/neum_Linux.Rds") {
    neum_Linux <- readRDS(file_names[i])
  }
  else if(file_names[i] == "./artifact/neum_macOS.Rds") {
    neum_macOS <- readRDS(file_names[i])
  }
  else if(file_names[i] == "./artifact/neum_Windows.Rds") {
    neum_Windows <- readRDS(file_names[i])
  }
}

cat(paste("Compare regular Linux with macOS:\n"))
print(all.equal(reg_Linux, reg_macOS, tolerance=1.0e-13))
cat(paste("\n\nCompare regular Linux with Windows:\n"))
print(all.equal(reg_Linux, reg_Windows, tolerance=1.0e-13))
cat(paste("\n\nCompare regular macOS with Windows:\n"))
print(all.equal(reg_macOS, reg_Windows, tolerance=1.0e-13))

cat(paste("\n\n\nCompare Neumaier Linux with macOS:\n"))
print(all.equal(neum_Linux, neum_macOS, tolerance=1.0e-13))
cat(paste("\n\nCompare Neumaier Linux with Windows:\n"))
print(all.equal(neum_Linux, neum_Windows, tolerance=1.0e-13))
cat(paste("\n\nCompare Neumaier macOS with Windows:\n"))
print(all.equal(neum_macOS, neum_Windows, tolerance=1.0e-13))