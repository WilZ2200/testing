###############################################################################
# 02_simulate_data.R
# Simulate realistic DCE response data for N=500 respondents
# Project: Health Preferences for Breast Cancer Screening
###############################################################################

library(MASS)
library(dplyr)
library(tidyr)
library(readr)

set.seed(42)

message("=== Simulating DCE Response Data ===")

# -------------------------------------------------------------------------
# Configuration
# -------------------------------------------------------------------------

N <- 500       # Number of respondents
n_sets <- 12   # Choice sets per respondent (1 block)
n_alts <- 3    # Alternatives per choice set (A, B, opt-out)

# -------------------------------------------------------------------------
# Load experimental design
# -------------------------------------------------------------------------

design <- read.csv("/Users/hengzhezhao/Desktop/testing/data/experimental_design.csv",
                   stringsAsFactors = TRUE)

# -------------------------------------------------------------------------
# True population parameters (conditional logit base)
# -------------------------------------------------------------------------

beta_true <- list(
  method_mri       =  0.45,
  method_us        = -0.30,
  freq_biennial    = -0.25,
  freq_triennial   = -0.55,
  cost             = -0.008,
  sensitivity      =  0.04,
  fpr              = -0.06,
  wait_1week       = -0.35,
  wait_3weeks      = -0.80,
  pain_mild        = -0.20,
  pain_moderate    = -0.55,
  optout           = -0.60
)

# -------------------------------------------------------------------------
# Heterogeneity parameters (SD for mixed logit)
# -------------------------------------------------------------------------

sigma_true <- list(
  method_mri    = 0.80,
  method_us     = 0.50,
  freq_biennial = 0.40,
  freq_triennial = 0.60,
  sensitivity   = 0.02,
  wait_3weeks   = 0.50,
  pain_moderate = 0.40
)

# -------------------------------------------------------------------------
# Latent class parameters (3 classes)
# -------------------------------------------------------------------------

class_probs <- c(0.35, 0.40, 0.25)
class_names <- c("Cost-conscious", "Accuracy-focused", "Convenience-seekers")

# Class-specific parameter adjustments (multipliers/overrides)
class_params <- list(
  # Class 1: Cost-conscious (35%)
  class1 = list(
    method_mri    =  0.20,
    method_us     = -0.15,
    freq_biennial = -0.30,
    freq_triennial = -0.60,
    cost          = -0.015,   # Higher cost sensitivity
    sensitivity   =  0.02,   # Lower sensitivity preference
    fpr           = -0.04,
    wait_1week    = -0.25,
    wait_3weeks   = -0.60,
    pain_mild     = -0.15,
    pain_moderate = -0.40,
    optout        = -0.40
  ),
  # Class 2: Accuracy-focused (40%)
  class2 = list(
    method_mri    =  0.70,
    method_us     = -0.10,
    freq_biennial = -0.15,
    freq_triennial = -0.40,
    cost          = -0.004,   # Lower cost sensitivity
    sensitivity   =  0.08,   # Higher sensitivity preference
    fpr           = -0.10,   # Higher FPR sensitivity
    wait_1week    = -0.20,
    wait_3weeks   = -0.50,
    pain_mild     = -0.10,
    pain_moderate = -0.35,
    optout        = -0.80
  ),
  # Class 3: Convenience-seekers (25%)
  class3 = list(
    method_mri    =  0.30,
    method_us     = -0.50,
    freq_biennial = -0.35,
    freq_triennial = -0.70,
    cost          = -0.006,
    sensitivity   =  0.03,
    fpr           = -0.05,
    wait_1week    = -0.60,   # Higher waiting time sensitivity
    wait_3weeks   = -1.20,   # Much higher
    pain_mild     = -0.40,   # Higher pain sensitivity
    pain_moderate = -0.90,   # Much higher
    optout        = -0.50
  )
)

# -------------------------------------------------------------------------
# Generate sociodemographic variables
# -------------------------------------------------------------------------

message("Generating sociodemographic variables...")

# Age: Normal(55, 8), truncated to 40-74
age_raw <- rnorm(N * 2, mean = 55, sd = 8)
age_raw <- age_raw[age_raw >= 40 & age_raw <= 74]
age <- round(age_raw[1:N])

# Income categories
income_probs <- c(0.15, 0.25, 0.35, 0.25)
income <- sample(c("<$30k", "$30-60k", "$60-100k", ">$100k"),
                 N, replace = TRUE, prob = income_probs)

# Education
edu_probs <- c(0.20, 0.30, 0.30, 0.20)
education <- sample(c("High school", "Some college", "Bachelor", "Graduate"),
                    N, replace = TRUE, prob = edu_probs)

# Family history of breast cancer
family_history <- rbinom(N, 1, 0.25)

# Prior screening experience
prior_screening <- rbinom(N, 1, 0.70)

# Insurance type
ins_probs <- c(0.55, 0.25, 0.10, 0.10)
insurance <- sample(c("Private", "Medicare", "Medicaid", "Uninsured"),
                    N, replace = TRUE, prob = ins_probs)

# Race/ethnicity
race_probs <- c(0.60, 0.15, 0.15, 0.07, 0.03)
race <- sample(c("White", "Black", "Hispanic", "Asian", "Other"),
               N, replace = TRUE, prob = race_probs)

# Create respondent dataframe
respondents <- data.frame(
  resp_id = 1:N,
  age = age,
  income = factor(income, levels = c("<$30k", "$30-60k", "$60-100k", ">$100k")),
  education = factor(education, levels = c("High school", "Some college",
                                            "Bachelor", "Graduate")),
  family_history = family_history,
  prior_screening = prior_screening,
  insurance = factor(insurance, levels = c("Private", "Medicare",
                                            "Medicaid", "Uninsured")),
  race = factor(race, levels = c("White", "Black", "Hispanic", "Asian", "Other")),
  stringsAsFactors = FALSE
)

# -------------------------------------------------------------------------
# Assign respondents to latent classes
# -------------------------------------------------------------------------

message("Assigning latent classes...")

# Class assignment based on class probabilities
# In practice, class membership may correlate with demographics
# Here we add slight demographic correlations
class_logit <- matrix(0, N, 3)
class_logit[, 1] <- 0  # Reference class

# Class 2 (Accuracy-focused): more likely if family history or higher education
class_logit[, 2] <- 0.3 + 0.5 * family_history +
  0.3 * (education == "Graduate") +
  0.2 * (education == "Bachelor")

# Class 3 (Convenience-seekers): more likely if younger or higher income
class_logit[, 3] <- -0.2 - 0.02 * (age - 55) +
  0.3 * (income == ">$100k") +
  0.2 * (income == "$60-100k")

# Convert to probabilities (softmax, scaled to match target proportions)
class_prob_individual <- exp(class_logit) / rowSums(exp(class_logit))

# Assign classes
latent_class <- numeric(N)
for (i in 1:N) {
  latent_class[i] <- sample(1:3, 1, prob = class_prob_individual[i, ])
}

respondents$latent_class <- latent_class
respondents$class_label <- class_names[latent_class]

message(paste("Class distribution:", paste(table(latent_class), collapse = ", ")))

# -------------------------------------------------------------------------
# Draw individual-level parameters (for mixed logit simulation)
# -------------------------------------------------------------------------

message("Drawing individual-level parameters...")

individual_params <- data.frame(resp_id = 1:N)

# For each respondent, start with class-specific means, add heterogeneity
for (i in 1:N) {
  cls <- latent_class[i]
  cls_pars <- class_params[[paste0("class", cls)]]

  # Draw random parameters with class-specific means and population SDs
  individual_params$method_mri[i]    <- rnorm(1, cls_pars$method_mri, sigma_true$method_mri)
  individual_params$method_us[i]     <- rnorm(1, cls_pars$method_us, sigma_true$method_us)
  individual_params$freq_biennial[i] <- rnorm(1, cls_pars$freq_biennial, sigma_true$freq_biennial)
  individual_params$freq_triennial[i] <- rnorm(1, cls_pars$freq_triennial, sigma_true$freq_triennial)
  individual_params$cost[i]          <- cls_pars$cost  # Cost kept fixed within class
  individual_params$sensitivity[i]   <- rnorm(1, cls_pars$sensitivity, sigma_true$sensitivity)
  individual_params$fpr[i]           <- cls_pars$fpr
  individual_params$wait_1week[i]    <- cls_pars$wait_1week
  individual_params$wait_3weeks[i]   <- rnorm(1, cls_pars$wait_3weeks, sigma_true$wait_3weeks)
  individual_params$pain_mild[i]     <- cls_pars$pain_mild
  individual_params$pain_moderate[i] <- rnorm(1, cls_pars$pain_moderate, sigma_true$pain_moderate)
  individual_params$optout[i]        <- cls_pars$optout
}

# -------------------------------------------------------------------------
# Generate choice data
# -------------------------------------------------------------------------

message("Generating choice data...")

# Get blocks and their choice sets
blocks <- sort(unique(design$block))
block_sets <- lapply(blocks, function(b) sort(unique(design$choice_set[design$block == b])))

# Assign each respondent to a block (roughly equal)
respondents$block <- rep(blocks, length.out = N)

# Helper: create design matrix row for a profile
create_x <- function(row) {
  x <- numeric(12)
  names(x) <- c("method_mri", "method_us", "freq_biennial", "freq_triennial",
                 "cost", "sensitivity", "fpr", "wait_1week", "wait_3weeks",
                 "pain_mild", "pain_moderate", "optout")

  # Check for opt-out
  if (row$alternative == "C") {
    x["optout"] <- 1
    return(x)
  }

  # Method dummies (ref = Mammography)
  x["method_mri"] <- as.numeric(row$method == "MRI")
  x["method_us"]  <- as.numeric(row$method == "Ultrasound")

  # Frequency dummies (ref = Annual)
  x["freq_biennial"]  <- as.numeric(row$frequency == "Biennial")
  x["freq_triennial"] <- as.numeric(row$frequency == "Triennial")

  # Continuous attributes
  x["cost"]        <- row$cost
  x["sensitivity"] <- row$sensitivity
  x["fpr"]         <- row$fpr

  # Waiting time dummies (ref = 1 day)
  x["wait_1week"]   <- as.numeric(row$waittime == "1 week")
  x["wait_3weeks"]  <- as.numeric(row$waittime == "3 weeks")

  # Pain dummies (ref = None)
  x["pain_mild"]     <- as.numeric(row$pain == "Mild")
  x["pain_moderate"] <- as.numeric(row$pain == "Moderate")

  # Not opt-out
  x["optout"] <- 0

  return(x)
}

# Storage for all choice observations
all_choices <- list()
obs_counter <- 0

for (i in 1:N) {
  resp_block <- respondents$block[i]
  choice_sets <- block_sets[[resp_block]]
  betas_i <- as.numeric(individual_params[i, -1])  # Drop resp_id

  for (cs in choice_sets) {
    cs_data <- design[design$choice_set == cs, ]

    # Calculate utility for each alternative
    utilities <- numeric(nrow(cs_data))
    for (j in 1:nrow(cs_data)) {
      x_j <- create_x(cs_data[j, ])
      v_j <- sum(betas_i * x_j)  # Deterministic utility
      # Add Gumbel error (Type I extreme value)
      e_j <- -log(-log(runif(1)))
      utilities[j] <- v_j + e_j
    }

    # Choose alternative with highest total utility
    chosen <- which.max(utilities)

    # Store each alternative as a row
    for (j in 1:nrow(cs_data)) {
      obs_counter <- obs_counter + 1
      all_choices[[obs_counter]] <- data.frame(
        resp_id = i,
        block = resp_block,
        choice_set = cs,
        alternative = cs_data$alternative[j],
        chosen = as.integer(j == chosen),
        method = as.character(cs_data$method[j]),
        frequency = as.character(cs_data$frequency[j]),
        cost = cs_data$cost[j],
        sensitivity = cs_data$sensitivity[j],
        fpr = cs_data$fpr[j],
        waittime = as.character(cs_data$waittime[j]),
        pain = as.character(cs_data$pain[j]),
        stringsAsFactors = FALSE
      )
    }
  }

  if (i %% 100 == 0) message(paste("  Processed respondent", i, "of", N))
}

# Combine all observations
dce_data <- bind_rows(all_choices)

message(paste("Total observations:", nrow(dce_data)))
message(paste("Respondents:", length(unique(dce_data$resp_id))))
message(paste("Choice sets per respondent:", n_sets))

# -------------------------------------------------------------------------
# Create dummy variables for analysis
# -------------------------------------------------------------------------

message("Creating analysis variables...")

dce_data <- dce_data %>%
  mutate(
    # Method dummies (ref = Mammography)
    method_mri = as.integer(method == "MRI"),
    method_us  = as.integer(method == "Ultrasound"),

    # Frequency dummies (ref = Annual)
    freq_biennial  = as.integer(frequency == "Biennial"),
    freq_triennial = as.integer(frequency == "Triennial"),

    # Waiting time dummies (ref = 1 day)
    wait_1week  = as.integer(waittime == "1 week"),
    wait_3weeks = as.integer(waittime == "3 weeks"),

    # Pain dummies (ref = None)
    pain_mild     = as.integer(pain == "Mild"),
    pain_moderate = as.integer(pain == "Moderate"),

    # Opt-out indicator
    optout = as.integer(alternative == "C"),

    # Create unique choice situation ID
    choice_id = paste0(resp_id, "_", choice_set)
  )

# -------------------------------------------------------------------------
# Merge with respondent demographics
# -------------------------------------------------------------------------

dce_full <- dce_data %>%
  left_join(respondents, by = c("resp_id", "block"))

# -------------------------------------------------------------------------
# Summary statistics
# -------------------------------------------------------------------------

message("\n=== Choice Summary ===")
choice_freq <- dce_data %>%
  filter(chosen == 1) %>%
  group_by(alternative) %>%
  summarise(n = n(), pct = round(n() / (N * n_sets) * 100, 1))
print(choice_freq)

message("\nChoice by method:")
method_freq <- dce_data %>%
  filter(chosen == 1, alternative != "C") %>%
  group_by(method) %>%
  summarise(n = n(), pct = round(n() / sum(dce_data$chosen[dce_data$alternative != "C"]) * 100, 1))
print(method_freq)

# -------------------------------------------------------------------------
# Save wide format (raw)
# -------------------------------------------------------------------------

write.csv(dce_full,
          "/Users/hengzhezhao/Desktop/testing/data/raw/simulated_dce_data.csv",
          row.names = FALSE)
message("\nRaw data saved to data/raw/simulated_dce_data.csv")

# -------------------------------------------------------------------------
# Create long format suitable for mlogit
# -------------------------------------------------------------------------

# mlogit requires a specific data structure:
# Each row = one alternative within a choice situation
# Variables: choice_id, alt, chosen, + attributes

dce_long <- dce_full %>%
  mutate(
    alt = case_when(
      alternative == "A" ~ 1L,
      alternative == "B" ~ 2L,
      alternative == "C" ~ 3L
    ),
    # Recode factors for analysis
    method_f = factor(method, levels = c("Mammography", "MRI", "Ultrasound",
                                          "No screening")),
    frequency_f = factor(frequency, levels = c("Annual", "Biennial",
                                                "Triennial", "None")),
    waittime_f = factor(waittime, levels = c("1 day", "1 week", "3 weeks", "None")),
    pain_f = factor(pain, levels = c("None", "Mild", "Moderate"))
  ) %>%
  arrange(resp_id, choice_set, alt)

write.csv(dce_long,
          "/Users/hengzhezhao/Desktop/testing/data/processed/dce_long.csv",
          row.names = FALSE)
message("Long-format data saved to data/processed/dce_long.csv")

# Save respondent data separately
write.csv(respondents,
          "/Users/hengzhezhao/Desktop/testing/data/processed/respondents.csv",
          row.names = FALSE)
message("Respondent data saved to data/processed/respondents.csv")

# Save individual parameters (ground truth) for validation
write.csv(individual_params,
          "/Users/hengzhezhao/Desktop/testing/data/processed/true_individual_params.csv",
          row.names = FALSE)

# Save true population parameters
saveRDS(list(beta_true = beta_true,
             sigma_true = sigma_true,
             class_params = class_params,
             class_probs = class_probs),
        "/Users/hengzhezhao/Desktop/testing/data/processed/true_parameters.rds")

message("\nData simulation complete.")
message(paste("Total rows in long format:", nrow(dce_long)))
message(paste("Respondents:", N))
message(paste("Choice observations per respondent:", n_sets))
