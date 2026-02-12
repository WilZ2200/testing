###############################################################################
# 01_experimental_design.R
# Generate D-optimal experimental design for the DCE
# Project: Health Preferences for Breast Cancer Screening
###############################################################################

library(AlgDesign)
set.seed(42)

message("=== Generating Experimental Design ===")

# -------------------------------------------------------------------------
# Define attribute levels
# -------------------------------------------------------------------------
# Attribute         | Levels                         | Coding
# method            | Mammography(ref), MRI, US      | Dummy
# frequency         | Annual(ref), Biennial, Tri     | Dummy
# cost              | 0, 50, 150, 300                | Continuous
# sensitivity       | 70, 85, 95                     | Continuous
# fpr               | 5, 10, 15                      | Continuous
# waittime          | 1day(ref), 1week, 3weeks       | Dummy
# pain              | None(ref), Mild, Moderate      | Dummy

# -------------------------------------------------------------------------
# Step 1: Generate full factorial design
# Full factorial: 3 x 3 x 4 x 3 x 3 x 3 x 3 = 2916 profiles
# -------------------------------------------------------------------------

full_fact <- gen.factorial(
  levels = c(3, 3, 4, 3, 3, 3, 3),
  varNames = c("method", "frequency", "cost", "sensitivity",
               "fpr", "waittime", "pain"),
  factors = "all"
)

message(paste("Full factorial profiles:", nrow(full_fact)))

# Recode to meaningful levels
full_fact$method <- factor(full_fact$method,
                           levels = 1:3,
                           labels = c("Mammography", "MRI", "Ultrasound"))

full_fact$frequency <- factor(full_fact$frequency,
                              levels = 1:3,
                              labels = c("Annual", "Biennial", "Triennial"))

cost_levels <- c(0, 50, 150, 300)
full_fact$cost <- cost_levels[as.numeric(full_fact$cost)]

sens_levels <- c(70, 85, 95)
full_fact$sensitivity <- sens_levels[as.numeric(full_fact$sensitivity)]

fpr_levels <- c(5, 10, 15)
full_fact$fpr <- fpr_levels[as.numeric(full_fact$fpr)]

full_fact$waittime <- factor(full_fact$waittime,
                             levels = 1:3,
                             labels = c("1 day", "1 week", "3 weeks"))

full_fact$pain <- factor(full_fact$pain,
                         levels = 1:3,
                         labels = c("None", "Mild", "Moderate"))

# -------------------------------------------------------------------------
# Step 2: Generate D-optimal fractional factorial design
# We need 72 profiles (36 choice sets x 2 alternatives)
# -------------------------------------------------------------------------

# For optFederov, we need to work with numeric coding
design_matrix <- full_fact
design_matrix$method_num <- as.numeric(design_matrix$method)
design_matrix$frequency_num <- as.numeric(design_matrix$frequency)
design_matrix$waittime_num <- as.numeric(design_matrix$waittime)
design_matrix$pain_num <- as.numeric(design_matrix$pain)

# For D-optimality with continuous attributes modeled as linear,
# the algorithm tends to select extreme levels. To ensure all levels appear
# (important for the DCE), treat all attributes as factors for design generation.
full_fact_design <- full_fact
full_fact_design$cost_f <- factor(full_fact_design$cost)
full_fact_design$sensitivity_f <- factor(full_fact_design$sensitivity)
full_fact_design$fpr_f <- factor(full_fact_design$fpr)

opt_design <- optFederov(
  ~ method + frequency + cost_f + sensitivity_f + fpr_f + waittime + pain,
  data = full_fact_design,
  nTrials = 72,
  criterion = "D",
  nRepeats = 50
)

message(paste("D-efficiency:", round(opt_design$Dea, 4)))

# Extract the selected profiles
design_profiles <- full_fact[opt_design$rows, ]
design_profiles$profile_id <- 1:72

# -------------------------------------------------------------------------
# Step 3: Pair profiles into choice sets
# -------------------------------------------------------------------------

# Create choice sets by pairing consecutive profiles
design_profiles$choice_set <- rep(1:36, each = 2)
design_profiles$alternative <- rep(c("A", "B"), times = 36)

# -------------------------------------------------------------------------
# Step 4: Assign to blocks
# 3 blocks of 12 choice sets each
# -------------------------------------------------------------------------

block_assignment <- rep(1:3, each = 12)
# Randomize block assignment
block_assignment <- sample(block_assignment)
design_profiles$block <- block_assignment[design_profiles$choice_set]

# -------------------------------------------------------------------------
# Step 5: Add opt-out alternative
# Each choice set has: Alternative A, Alternative B, Opt-out (no screening)
# -------------------------------------------------------------------------

all_method_levels <- c("Mammography", "MRI", "Ultrasound", "No screening")
all_freq_levels <- c("Annual", "Biennial", "Triennial", "None")
all_wait_levels <- c("1 day", "1 week", "3 weeks", "None")
all_pain_levels <- c("None", "Mild", "Moderate")  # "None" already exists

optout_rows <- data.frame(
  method = factor(rep("No screening", 36), levels = all_method_levels),
  frequency = factor(rep("None", 36), levels = all_freq_levels),
  cost = 0,
  sensitivity = 0,
  fpr = 0,
  waittime = factor(rep("None", 36), levels = all_wait_levels),
  pain = factor(rep("None", 36), levels = all_pain_levels),
  profile_id = 73:108,
  choice_set = 1:36,
  alternative = "C",
  block = block_assignment
)

# For the final design, update factor levels to include opt-out
design_profiles$method <- factor(design_profiles$method, levels = all_method_levels)
design_profiles$frequency <- factor(design_profiles$frequency, levels = all_freq_levels)
design_profiles$waittime <- factor(design_profiles$waittime, levels = all_wait_levels)
design_profiles$pain <- factor(design_profiles$pain, levels = all_pain_levels)

design_final <- rbind(design_profiles, optout_rows)
design_final <- design_final[order(design_final$choice_set, design_final$alternative), ]
rownames(design_final) <- NULL

# -------------------------------------------------------------------------
# Step 6: Calculate and report design diagnostics
# -------------------------------------------------------------------------

message("\n=== Design Summary ===")
message(paste("Total profiles:", nrow(design_profiles)))
message(paste("Choice sets:", max(design_profiles$choice_set)))
message(paste("Alternatives per set: 2 + opt-out"))
message(paste("Blocks:", max(design_profiles$block)))
message(paste("Choice sets per block:", table(block_assignment)[1]))
message(paste("D-efficiency:", round(opt_design$Dea, 4)))

# Level balance check
message("\n--- Level Balance ---")
for (attr_name in c("method", "frequency", "cost", "sensitivity",
                     "fpr", "waittime", "pain")) {
  message(paste("\n", attr_name, ":"))
  print(table(design_profiles[[attr_name]]))
}

# -------------------------------------------------------------------------
# Step 7: Save the design
# -------------------------------------------------------------------------

write.csv(design_final,
          "/Users/hengzhezhao/Desktop/testing/data/experimental_design.csv",
          row.names = FALSE)

message("\nExperimental design saved to data/experimental_design.csv")

# Also save design metadata
design_meta <- list(
  n_profiles = nrow(design_profiles),
  n_choice_sets = 36,
  n_alternatives = 3,
  n_blocks = 3,
  d_efficiency = opt_design$Dea,
  seed = 42
)

saveRDS(design_meta,
        "/Users/hengzhezhao/Desktop/testing/data/processed/design_metadata.rds")

message("Design generation complete.")
