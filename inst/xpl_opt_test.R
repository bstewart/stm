#!/usr/bin/env Rscript
library(stringr)

# Get all the file names
location <- "./artifact/"
#location <- "~/sociology/stm/git/sumcpp/stm/inst/docopt/"
files <- list.files(path = location)

# Variables
num_files <- length(files)
os_names <- character()
bad_docs <- numeric()

# Loop through all the file names that have naming convention: 
# methodType_tolerance_osType.Rds
# and extract different operating systems
for(i in 1:num_files) {
  fnmx <- str_split_fixed(files[i], "_|\\.", 4)
  os_names <- cbind(os_names, fnmx[1,3])
}

# Because the method, and file extension is the same within a run we can do this
m_type <- fnmx[1,1]
tolerance <- fnmx[1,2]
f_ext <- fnmx[1,4]

# Remove duplicates
os_names <- unique(os_names, MARGIN = 2)

# Generate all pair combinations of os-s.
comb_os <- combn(os_names, 2)

# Check equality for all combinations of os-s for each topic
for(j in 1:ncol(comb_os)) {
  cat(paste("Comparing", comb_os[1,j], "vs", comb_os[2,j], "at", tolerance, "\n"))
  file_name1 <- str_c(location, m_type, "_", tolerance, "_", comb_os[1,j], ".", f_ext)
  file_name2 <- str_c(location, m_type, "_", tolerance, "_", comb_os[2,j], ".", f_ext)
  os1 <- readRDS(file_name1)
  os2 <- readRDS(file_name2)
  
  # Search for bad document
  for(i in 1:nrow(os1)) {
    if(!isTRUE(all.equal(os1$sol_norm2[i], os2$sol_norm2[i]))) {
      bad_docs <- cbind(bad_docs, i)
    }
  }
  
  cat(paste("Number of bad docs:", length(bad_docs), "for method:", m_type, "at", tolerance, ".\n"))
  #print(bad_docs)
}
