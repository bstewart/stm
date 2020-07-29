#!/usr/bin/env Rscript

# Get command line argument for operating system
args = commandArgs(trailingOnly=TRUE)

file_names = c("./artifact/reg_Linux.Rds", "./artifact/reg_macOS.Rds", "./artifact/reg_Windows.Rds", "./artifact/neum_Linux.Rds", "./artifact/neum_macOS.Rds", "./artifact/neum_Windows.Rds")

tol = 1.0e-9

for(i in 1:length(file_names)) {
  if(file_names[i] == "./artifact/reg_Linux.Rds") {
    cat(paste("Reading", file_names[i], "...\n"))
    reg_Linux <- readRDS(file_names[i])
  }
  else if(file_names[i] == "./artifact/reg_macOS.Rds") {
    cat(paste("Reading", file_names[i], "...\n"))
    reg_macOS <- readRDS(file_names[i])
  }
  else if(file_names[i] == "./artifact/reg_Windows.Rds") {
    cat(paste("Reading", file_names[i], "...\n"))
    reg_Windows <- readRDS(file_names[i])
  }
  else if(file_names[i] == "./artifact/neum_Linux.Rds") {
    cat(paste("Reading", file_names[i], "...\n"))
    neum_Linux <- readRDS(file_names[i])
  }
  else if(file_names[i] == "./artifact/neum_macOS.Rds") {
    cat(paste("Reading", file_names[i], "...\n"))
    neum_macOS <- readRDS(file_names[i])
  }
  else if(file_names[i] == "./artifact/neum_Windows.Rds") {
    cat(paste("Reading", file_names[i], "...\n"))
    neum_Windows <- readRDS(file_names[i])
  }
}

cat(paste("Compare regular Linux with macOS:\n"))
print(all.equal(reg_Linux, reg_macOS, tolerance=tol))
cat(paste("\n\nCompare regular Linux with Windows:\n"))
print(all.equal(reg_Linux, reg_Windows, tolerance=tol))
cat(paste("\n\nCompare regular macOS with Windows:\n"))
print(all.equal(reg_macOS, reg_Windows, tolerance=tol))

cat(paste("\n\n\nCompare Neumaier Linux with macOS:\n"))
print(all.equal(neum_Linux, neum_macOS, tolerance=tol))
cat(paste("\n\nCompare Neumaier Linux with Windows:\n"))
print(all.equal(neum_Linux, neum_Windows, tolerance=tol))
cat(paste("\n\nCompare Neumaier macOS with Windows:\n"))
print(all.equal(neum_macOS, neum_Windows, tolerance=tol))

cat(paste("\n\n\nCompare Neumaier Linux with regular Linux:\n"))
print(all.equal(neum_Linux, reg_Linux, tolerance=tol))
cat(paste("\n\nCompare Neumaier Windows with regular Windows:\n"))
print(all.equal(neum_Windows, reg_Windows, tolerance=tol))
cat(paste("\n\nCompare Neumaier macOS with regular macOS:\n"))
print(all.equal(neum_macOS, reg_macOS, tolerance=tol))