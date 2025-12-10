# Complete Setup and Run Script
# Run this to set everything up and compare SVI vs EM

cat("=======================================================\n")
cat("STM with SVI - Complete Setup and Comparison\n")
cat("=======================================================\n\n")

# Step 1: Update NAMESPACE
cat("Step 1: Updating NAMESPACE...\n")
if(!requireNamespace("roxygen2", quietly=TRUE)) {
  install.packages("roxygen2", repos="https://cloud.r-project.org")
}
roxygen2::roxygenise()

# Verify
ns <- readLines("NAMESPACE")
if(any(grepl("export\\(stm_svi\\)", ns))) {
  cat("✓ stm_svi exported in NAMESPACE\n\n")
} else {
  cat("✗ Warning: stm_svi not in NAMESPACE\n\n")
}

# Step 2: Instructions for user
cat("=======================================================\n")
cat("IMPORTANT: You must now RESTART R!\n")
cat("=======================================================\n\n")

cat("Please do the following:\n")
cat("1. Restart R session:\n")
cat("   - In RStudio: Session -> Restart R (Ctrl+Shift+F10)\n")
cat("   - Or run: .rs.restartR()\n\n")

cat("2. After restarting, run:\n")
cat("   library(stm)\n")
cat("   source('quick_compare_installed.R')\n\n")

cat("=======================================================\n")
cat("Alternative: Use this two-step process\n")
cat("=======================================================\n\n")

cat("Step A (this file):\n")
cat("  source('SETUP_AND_RUN.R')  # Updates NAMESPACE\n\n")

cat("Step B (after restart):\n")
cat("  .rs.restartR()  # Restart R\n")
cat("  library(stm)\n")
cat("  source('quick_compare_installed.R')  # Run comparison\n\n")

cat("Or use the automated version (if you have devtools):\n")
cat("  source('build_package.R')  # Does everything including install\n\n")

cat("=======================================================\n")
