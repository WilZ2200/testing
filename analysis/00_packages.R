###############################################################################
# 00_packages.R
# Install and load all required packages for DCE analysis
# Project: Health Preferences for Breast Cancer Screening
###############################################################################

required_packages <- c(
  "mlogit",        # Multinomial logit models
  "gmnl",          # Generalized multinomial logit (mixed logit and latent class)
  "AlgDesign",     # Experimental design generation
  "survival",      # Conditional logit via clogit
  "MASS",          # Multivariate normal simulation
  "ggplot2",       # Plotting
  "dplyr",         # Data manipulation
  "tidyr",         # Data reshaping
  "stargazer",     # LaTeX table output
  "xtable",        # Additional LaTeX tables
  "scales",        # Scale functions for plots
  "RColorBrewer",  # Color palettes
  "gridExtra",     # Multiple plots
  "broom",         # Tidy model output
  "patchwork",     # Combine ggplots
  "haven",         # Read/write data formats
  "readr",         # CSV reading
  "stringr",       # String manipulation
  "forcats",       # Factor manipulation
  "kableExtra",    # Better tables
  "modelsummary"   # Model comparison tables
)

# Install missing packages
install_if_missing <- function(pkg) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    message(paste("Installing package:", pkg))
    install.packages(pkg, repos = "https://cran.r-project.org", quiet = TRUE)
    library(pkg, character.only = TRUE)
  }
}

invisible(sapply(required_packages, install_if_missing))

message("All packages loaded successfully.")
