###############################################################################
# run_all.R
# Master script: run all analysis scripts in order
# Project: Health Preferences for Breast Cancer Screening
###############################################################################

# Record start time
start_time <- Sys.time()

message("============================================================")
message("  DCE Analysis Pipeline")
message("  Health Preferences for Breast Cancer Screening")
message(paste("  Started:", start_time))
message("============================================================")

# Set working directory to analysis folder
setwd("/Users/hengzhezhao/Desktop/testing/analysis")

# Set global seed
set.seed(42)

# Step 0: Install and load packages
message("\n>>> Step 0: Loading packages...")
source("00_packages.R")

# Step 1: Generate experimental design
message("\n>>> Step 1: Generating experimental design...")
source("01_experimental_design.R")

# Step 2: Simulate data
message("\n>>> Step 2: Simulating DCE data...")
source("02_simulate_data.R")

# Step 3: Run analysis models
message("\n>>> Step 3: Running analysis models...")
source("03_analysis.R")

# Step 4: Generate figures
message("\n>>> Step 4: Generating figures...")
source("04_figures.R")

# Step 5: Generate tables
message("\n>>> Step 5: Generating LaTeX tables...")
source("05_tables.R")

# Report completion
end_time <- Sys.time()
elapsed <- difftime(end_time, start_time, units = "mins")

message("\n============================================================")
message("  Pipeline Complete")
message(paste("  Ended:", end_time))
message(paste("  Total time:", round(as.numeric(elapsed), 1), "minutes"))
message("============================================================")

# List output files
message("\nOutput files:")
message("  Data:")
message("    data/experimental_design.csv")
message("    data/raw/simulated_dce_data.csv")
message("    data/processed/dce_long.csv")
message("    data/processed/respondents.csv")
message("  Figures:")
for (f in list.files("/Users/hengzhezhao/Desktop/testing/figures/")) {
  message(paste("   ", f))
}
message("  Tables:")
for (f in list.files("/Users/hengzhezhao/Desktop/testing/tables/")) {
  message(paste("   ", f))
}
