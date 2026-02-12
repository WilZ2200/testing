###############################################################################
# 04_figures.R
# Generate publication-quality figures for DCE analysis
# Project: Health Preferences for Breast Cancer Screening
###############################################################################

library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(RColorBrewer)
library(scales)
library(gridExtra)
library(MASS)

# Ensure dplyr::select is used (not MASS::select)
select <- dplyr::select

set.seed(42)

message("=== Generating Figures ===")

# -------------------------------------------------------------------------
# Load results
# -------------------------------------------------------------------------

cl_model <- readRDS("/Users/hengzhezhao/Desktop/testing/data/processed/cl_model.rds")
cl_coefs <- readRDS("/Users/hengzhezhao/Desktop/testing/data/processed/cl_coefs.rds")
mxl_model <- readRDS("/Users/hengzhezhao/Desktop/testing/data/processed/mxl_model.rds")
mxl_coefs <- readRDS("/Users/hengzhezhao/Desktop/testing/data/processed/mxl_coefs.rds")
lc_best <- readRDS("/Users/hengzhezhao/Desktop/testing/data/processed/lc_best.rds")
wtp_results <- readRDS("/Users/hengzhezhao/Desktop/testing/data/processed/wtp_results.rds")
dce_long <- read.csv("/Users/hengzhezhao/Desktop/testing/data/processed/dce_long.csv",
                     stringsAsFactors = TRUE)

# -------------------------------------------------------------------------
# Theme setup
# -------------------------------------------------------------------------

theme_dce <- theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_line(color = "grey90"),
    panel.grid.major.x = element_blank(),
    axis.text = element_text(color = "black"),
    axis.title = element_text(face = "bold"),
    plot.title = element_text(face = "bold", hjust = 0),
    plot.subtitle = element_text(color = "grey40"),
    legend.position = "bottom",
    strip.text = element_text(face = "bold")
  )

# Colorblind-friendly palette
cb_palette <- c("#0072B2", "#D55E00", "#009E73", "#CC79A7",
                "#F0E442", "#56B4E9", "#E69F00")

# =========================================================================
# Figure 1: Coefficient Plot (Forest Plot)
# =========================================================================

message("Generating Figure 1: Coefficient plot...")

# Prepare CL data
cl_plot_data <- cl_coefs %>%
  filter(variable != "optout") %>%
  mutate(
    label = c("MRI vs. Mammography", "Ultrasound vs. Mammography",
              "Biennial vs. Annual", "Triennial vs. Annual",
              "Cost (per $1)", "Sensitivity (per 1pp)",
              "False positive rate (per 1pp)",
              "1 week vs. 1 day wait", "3 weeks vs. 1 day wait",
              "Mild vs. No pain", "Moderate vs. No pain"),
    attribute_group = c("Method", "Method",
                        "Frequency", "Frequency",
                        "Cost", "Sensitivity", "FPR",
                        "Wait time", "Wait time",
                        "Pain", "Pain"),
    model = "Conditional Logit"
  )

# Prepare MXL means
mxl_means <- mxl_coefs %>%
  filter(!grepl("^sd\\.", variable)) %>%
  filter(variable != "optout" & variable != "cost" & variable != "fpr") %>%
  mutate(ci_lower = estimate - 1.96 * se,
         ci_upper = estimate + 1.96 * se)

# Match labels for MXL
mxl_label_map <- c(
  "method_mri" = "MRI vs. Mammography",
  "method_us" = "Ultrasound vs. Mammography",
  "freq_biennial" = "Biennial vs. Annual",
  "freq_triennial" = "Triennial vs. Annual",
  "sensitivity" = "Sensitivity (per 1pp)",
  "wait_1week" = "1 week vs. 1 day wait",
  "wait_3weeks" = "3 weeks vs. 1 day wait",
  "pain_mild" = "Mild vs. No pain",
  "pain_moderate" = "Moderate vs. No pain"
)

mxl_plot_data <- mxl_means %>%
  filter(variable %in% names(mxl_label_map)) %>%
  mutate(
    label = mxl_label_map[variable],
    model = "Mixed Logit"
  )

# Also get CL data for matching variables
cl_match <- cl_plot_data %>%
  filter(label %in% mxl_plot_data$label)

# Combine
coef_plot_data <- bind_rows(
  cl_match %>% select(label, estimate, ci_lower, ci_upper, model, attribute_group),
  mxl_plot_data %>%
    mutate(attribute_group = case_when(
      grepl("MRI|Ultrasound", label) ~ "Method",
      grepl("ennial", label) ~ "Frequency",
      grepl("Sensitivity", label) ~ "Sensitivity",
      grepl("week|wait", label) ~ "Wait time",
      grepl("pain|Pain", label) ~ "Pain",
      TRUE ~ "Other"
    )) %>%
    select(label, estimate, ci_lower, ci_upper, model, attribute_group)
)

# Order labels
label_order <- rev(c("MRI vs. Mammography", "Ultrasound vs. Mammography",
                      "Biennial vs. Annual", "Triennial vs. Annual",
                      "Sensitivity (per 1pp)",
                      "1 week vs. 1 day wait", "3 weeks vs. 1 day wait",
                      "Mild vs. No pain", "Moderate vs. No pain"))

coef_plot_data$label <- factor(coef_plot_data$label, levels = label_order)

p1 <- ggplot(coef_plot_data, aes(x = estimate, y = label, color = model, shape = model)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
  geom_point(position = position_dodge(width = 0.5), size = 3) +
  geom_errorbarh(aes(xmin = ci_lower, xmax = ci_upper),
                 position = position_dodge(width = 0.5), height = 0.2) +
  scale_color_manual(values = c("Conditional Logit" = cb_palette[1],
                                "Mixed Logit" = cb_palette[2])) +
  scale_shape_manual(values = c("Conditional Logit" = 16, "Mixed Logit" = 17)) +
  labs(
    title = "Estimated Coefficients from DCE Models",
    subtitle = "Point estimates with 95% confidence intervals",
    x = "Coefficient Estimate",
    y = "",
    color = "Model",
    shape = "Model"
  ) +
  theme_dce +
  theme(legend.position = "bottom")

ggsave("/Users/hengzhezhao/Desktop/testing/figures/coefficient_plot.pdf",
       p1, width = 10, height = 7)
ggsave("/Users/hengzhezhao/Desktop/testing/figures/coefficient_plot.png",
       p1, width = 10, height = 7, dpi = 300)

message("  Saved coefficient_plot.pdf and .png")

# =========================================================================
# Figure 2: WTP Plot
# =========================================================================

message("Generating Figure 2: WTP plot...")

wtp_kr <- wtp_results$kr

# Order by WTP magnitude
wtp_kr$attribute <- factor(wtp_kr$attribute,
                           levels = wtp_kr$attribute[order(wtp_kr$wtp_mean)])

p2 <- ggplot(wtp_kr, aes(x = wtp_mean, y = attribute)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
  geom_point(size = 3, color = cb_palette[1]) +
  geom_errorbarh(aes(xmin = ci_lower, xmax = ci_upper),
                 height = 0.3, color = cb_palette[1]) +
  labs(
    title = "Willingness-to-Pay Estimates",
    subtitle = "Krinsky-Robb 95% confidence intervals (CL model)",
    x = "WTP (USD)",
    y = ""
  ) +
  theme_dce +
  scale_x_continuous(labels = dollar_format())

ggsave("/Users/hengzhezhao/Desktop/testing/figures/wtp_plot.pdf",
       p2, width = 10, height = 7)
ggsave("/Users/hengzhezhao/Desktop/testing/figures/wtp_plot.png",
       p2, width = 10, height = 7, dpi = 300)

message("  Saved wtp_plot.pdf and .png")

# =========================================================================
# Figure 3: Latent Class Profile Plot
# =========================================================================

message("Generating Figure 3: Latent class profiles...")

# LC results use custom format: lc_best$betas (matrix), lc_best$pi_q, lc_best$Q
if (!is.null(lc_best) && !is.null(lc_best$betas)) {
  Q <- lc_best$Q

  param_names <- rownames(lc_best$betas)
  pretty_names <- c("MRI", "Ultrasound", "Biennial", "Triennial",
                    "Cost", "Sensitivity", "FPR", "Wait 1wk", "Wait 3wks",
                    "Pain mild", "Pain moderate", "Opt-out")

  lc_class_coefs <- data.frame()
  for (q in 1:Q) {
    for (k in seq_along(param_names)) {
      lc_class_coefs <- rbind(lc_class_coefs, data.frame(
        class = paste0("Class ", q, " (", round(lc_best$pi_q[q] * 100, 0), "%)"),
        parameter = pretty_names[k],
        estimate = lc_best$betas[param_names[k], q],
        stringsAsFactors = FALSE
      ))
    }
  }

  lc_class_coefs$parameter <- factor(lc_class_coefs$parameter, levels = pretty_names)

  p3 <- ggplot(lc_class_coefs %>% filter(parameter != "Opt-out"),
               aes(x = parameter, y = estimate, fill = class)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.8),
             width = 0.7) +
    geom_hline(yintercept = 0, linetype = "solid", color = "grey40") +
    scale_fill_manual(values = cb_palette[1:Q]) +
    labs(
      title = "Latent Class Coefficient Profiles",
      subtitle = paste0(Q, "-class model"),
      x = "",
      y = "Coefficient Estimate",
      fill = "Latent Class"
    ) +
    theme_dce +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  ggsave("/Users/hengzhezhao/Desktop/testing/figures/latent_class_profiles.pdf",
         p3, width = 12, height = 7)
  ggsave("/Users/hengzhezhao/Desktop/testing/figures/latent_class_profiles.png",
         p3, width = 12, height = 7, dpi = 300)

  message("  Saved latent_class_profiles.pdf and .png")
} else {
  message("  Warning: No LC results available; skipping LC profile plot")
}

# =========================================================================
# Figure 4: Predicted Probability Plot
# =========================================================================

message("Generating Figure 4: Predicted probabilities...")

# Use CL model to predict screening uptake across cost and sensitivity scenarios
cl_beta <- coef(cl_model)

# Create scenarios: vary cost and sensitivity, fix other attributes at reference
cost_range <- seq(0, 300, by = 10)
sens_levels <- c(70, 85, 95)

pred_data <- expand.grid(cost = cost_range, sensitivity = sens_levels)
pred_data$prob <- NA

for (i in 1:nrow(pred_data)) {
  # Utility for screening alternative (Mammography, Annual, reference waittime/pain)
  v_screen <- cl_beta["cost"] * pred_data$cost[i] +
    cl_beta["sensitivity"] * pred_data$sensitivity[i] +
    cl_beta["fpr"] * 10  # Fix FPR at 10%

  # Utility for opt-out
  v_optout <- cl_beta["optout"]

  # Probability of choosing screening (logit)
  # In binary choice: P = exp(V_screen) / (exp(V_screen) + exp(V_optout))
  pred_data$prob[i] <- exp(v_screen) / (exp(v_screen) + exp(v_optout))
}

pred_data$sensitivity_label <- paste0(pred_data$sensitivity, "% sensitivity")

p4 <- ggplot(pred_data, aes(x = cost, y = prob, color = sensitivity_label)) +
  geom_line(linewidth = 1.2) +
  scale_color_manual(values = cb_palette[1:3]) +
  scale_y_continuous(labels = percent_format(), limits = c(0, 1)) +
  scale_x_continuous(labels = dollar_format()) +
  labs(
    title = "Predicted Screening Uptake by Cost and Sensitivity",
    subtitle = "Mammography, annual, reference levels for other attributes",
    x = "Out-of-Pocket Cost",
    y = "Predicted Probability of Screening Uptake",
    color = ""
  ) +
  theme_dce

ggsave("/Users/hengzhezhao/Desktop/testing/figures/predicted_probabilities.pdf",
       p4, width = 10, height = 7)
ggsave("/Users/hengzhezhao/Desktop/testing/figures/predicted_probabilities.png",
       p4, width = 10, height = 7, dpi = 300)

message("  Saved predicted_probabilities.pdf and .png")

# =========================================================================
# Figure 5: Heterogeneity Plot (MXL Random Parameter Distributions)
# =========================================================================

message("Generating Figure 5: Heterogeneity plot...")

# Extract mean and SD of random parameters from MXL
mxl_all <- coef(mxl_model)
mxl_names <- names(mxl_all)

# Random parameter names
rpar_names <- c("method_mri", "method_us", "freq_biennial", "freq_triennial",
                "sensitivity", "wait_1week", "wait_3weeks",
                "pain_mild", "pain_moderate")

rpar_labels <- c("MRI", "Ultrasound", "Biennial", "Triennial",
                 "Sensitivity", "Wait 1 week", "Wait 3 weeks",
                 "Pain mild", "Pain moderate")

# Build distributions
het_data <- data.frame()
x_range <- seq(-3, 3, length.out = 200)

for (k in seq_along(rpar_names)) {
  pname <- rpar_names[k]
  sd_name <- paste0("sd.", pname)

  if (pname %in% mxl_names && sd_name %in% mxl_names) {
    mu_k <- mxl_all[pname]
    sd_k <- abs(mxl_all[sd_name])

    x_vals <- mu_k + x_range * sd_k
    y_vals <- dnorm(x_vals, mean = mu_k, sd = sd_k)

    het_data <- rbind(het_data, data.frame(
      parameter = rpar_labels[k],
      x = x_vals,
      density = y_vals,
      stringsAsFactors = FALSE
    ))
  }
}

if (nrow(het_data) > 0) {
  p5 <- ggplot(het_data, aes(x = x, y = density)) +
    geom_line(color = cb_palette[1], linewidth = 0.8) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
    geom_area(alpha = 0.2, fill = cb_palette[1]) +
    facet_wrap(~ parameter, scales = "free", ncol = 3) +
    labs(
      title = "Distribution of Random Parameters (Mixed Logit)",
      subtitle = "Normal distributions showing preference heterogeneity",
      x = "Parameter Value",
      y = "Density"
    ) +
    theme_dce

  ggsave("/Users/hengzhezhao/Desktop/testing/figures/heterogeneity_plot.pdf",
         p5, width = 12, height = 10)
  ggsave("/Users/hengzhezhao/Desktop/testing/figures/heterogeneity_plot.png",
         p5, width = 12, height = 10, dpi = 300)

  message("  Saved heterogeneity_plot.pdf and .png")
} else {
  message("  Warning: Could not extract MXL random parameters; skipping heterogeneity plot")
}

message("\n=== All figures generated ===")
