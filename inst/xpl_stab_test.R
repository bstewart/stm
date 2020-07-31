#!/usr/bin/env Rscript
library(stringr)

# Get all the file names
location <- "./artifact/"
files <- list.files(path = location)
print(files)

# Variables
num_files <- length(files)
os_names <- character()
num_topics <- character()

# Loop through all the file names that have naming convention: 
# numTopics_sumType_methodType_osType.Rds
# and extract number of topics and different operating systems
for(i in 1:num_files) {
  fnmx <- str_split_fixed(files[i], "_|\\.", 5)
  print(fnmx)
  num_topics <- cbind(num_topics, fnmx[1,1])
  os_names <- cbind(os_names, fnmx[1,4])
}

# Because the summation type, the method, and file extension is the same within a run we can do this
s_type <- fnmx[1,2]
m_type <- fnmx[1,3]
f_ext <- fnmx[1,5]

# Remove duplicates
num_topics <- unique(num_topics, MARGIN = 2)
os_names <- unique(os_names, MARGIN = 2)

# Generate all pair combinations of os-s.
comb_os <- combn(os_names, 2)

# Check equality for all topics and combinations of os-s for each topic
for(i in 1:length(num_topics)) {
  cat(paste("\n\nChecking equality for K =", num_topics[i], "\n"))
  for(j in 1:ncol(comb_os)) {
    cat(paste("Comparing", comb_os[1,j], "vs", comb_os[2,j], "\n"))
    file_name1 <- str_c(location, num_topics[i], "_", s_type, "_", m_type, "_", comb_os[1,j], ".", f_ext)
    file_name2 <- str_c(location, num_topics[i], "_", s_type, "_", m_type, "_", comb_os[2,j], ".", f_ext)
    os1 <- readRDS(file_name1)
    os2 <- readRDS(file_name2)
    comp <- all.equal(os1, os2)
    if(length(comp) < 8) {
      cat("SUCCEEDED!!!\n")
    }
    else {
      cat("FAILED!!!\n")
    }
    print(comp)
  }
}