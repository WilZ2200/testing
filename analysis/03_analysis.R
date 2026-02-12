###############################################################################
# 03_analysis.R
# Full DCE analysis pipeline: CL, MXL, LC models, WTP, subgroup analysis
# Project: Health Preferences for Breast Cancer Screening
###############################################################################

library(mlogit)
library(gmnl)
library(dplyr)
library(tidyr)
library(MASS)
library(readr)
library(broom)

set.seed(42)

message("=== Starting DCE Analysis ===")

# -------------------------------------------------------------------------
# Load data
# -------------------------------------------------------------------------

dce_long <- read.csv("/Users/hengzhezhao/Desktop/testing/data/processed/dce_long.csv",
                     stringsAsFactors = TRUE)
respondents <- read.csv("/Users/hengzhezhao/Desktop/testing/data/processed/respondents.csv",
                        stringsAsFactors = TRUE)

message(paste("Loaded", nrow(dce_long), "observations from",
              length(unique(dce_long$resp_id)), "respondents"))

# =========================================================================
# 3.1 Descriptive Statistics
# =========================================================================

message("\n=== 3.1 Descriptive Statistics ===")

# Sample characteristics
message("\n--- Sample Characteristics ---")
message(paste("N =", length(unique(respondents$resp_id))))
message(paste("Age: mean =", round(mean(respondents$age), 1),
              ", SD =", round(sd(respondents$age), 1)))

# Demographics tables
demo_tables <- list(
  age_summary = respondents %>%
    summarise(mean = mean(age), sd = sd(age), min = min(age), max = max(age)),
  income = respondents %>% count(income) %>% mutate(pct = round(n/sum(n)*100, 1)),
  education = respondents %>% count(education) %>% mutate(pct = round(n/sum(n)*100, 1)),
  family_history = respondents %>% count(family_history) %>% mutate(pct = round(n/sum(n)*100, 1)),
  prior_screening = respondents %>% count(prior_screening) %>% mutate(pct = round(n/sum(n)*100, 1)),
  insurance = respondents %>% count(insurance) %>% mutate(pct = round(n/sum(n)*100, 1)),
  race = respondents %>% count(race) %>% mutate(pct = round(n/sum(n)*100, 1))
)

# Choice frequencies
choice_summary <- dce_long %>%
  filter(chosen == 1) %>%
  group_by(alternative) %>%
  summarise(n = n()) %>%
  mutate(pct = round(n / sum(n) * 100, 1))

message("\nChoice frequencies:")
print(choice_summary)

# Attribute level frequencies in chosen alternatives
attr_freq <- dce_long %>%
  filter(chosen == 1, alternative != "C") %>%
  group_by(method) %>%
  summarise(n = n()) %>%
  mutate(pct = round(n / sum(n) * 100, 1))

message("\nChosen screening method:")
print(attr_freq)

# Save descriptive results
saveRDS(demo_tables,
        "/Users/hengzhezhao/Desktop/testing/data/processed/demo_tables.rds")
saveRDS(choice_summary,
        "/Users/hengzhezhao/Desktop/testing/data/processed/choice_summary.rds")

# =========================================================================
# 3.2 Conditional Logit Model
# =========================================================================

message("\n=== 3.2 Conditional Logit Model ===")

# Prepare mlogit data
# Create proper choice index
dce_long$chid <- dce_long$choice_id
dce_long$alt_id <- dce_long$alt

mlogit_data <- dfidx(dce_long,
                     choice = "chosen",
                     idx = list(c("chid", "resp_id"), "alt_id"),
                     drop.index = FALSE)

# Estimate conditional logit
# Using generic attributes (not alternative-specific)
cl_formula <- chosen ~ method_mri + method_us + freq_biennial + freq_triennial +
  cost + sensitivity + fpr + wait_1week + wait_3weeks +
  pain_mild + pain_moderate + optout | 0

cl_model <- mlogit(cl_formula, data = mlogit_data)

message("\n--- Conditional Logit Results ---")
print(summary(cl_model))

# Extract results
cl_coefs <- data.frame(
  variable = names(coef(cl_model)),
  estimate = coef(cl_model),
  se = sqrt(diag(vcov(cl_model))),
  stringsAsFactors = FALSE
)
cl_coefs$z_value <- cl_coefs$estimate / cl_coefs$se
cl_coefs$p_value <- 2 * pnorm(-abs(cl_coefs$z_value))
cl_coefs$ci_lower <- cl_coefs$estimate - 1.96 * cl_coefs$se
cl_coefs$ci_upper <- cl_coefs$estimate + 1.96 * cl_coefs$se
cl_coefs$sig <- ifelse(cl_coefs$p_value < 0.001, "***",
                ifelse(cl_coefs$p_value < 0.01, "**",
                ifelse(cl_coefs$p_value < 0.05, "*",
                ifelse(cl_coefs$p_value < 0.1, ".", ""))))
rownames(cl_coefs) <- NULL

message("\nCL Coefficient Table:")
print(cl_coefs)

# Model fit
cl_ll <- logLik(cl_model)
n_choice_obs <- nrow(dce_long) / 3  # Number of choice situations
cl_ll_null <- -log(3) * n_choice_obs  # Null model: equal probability for 3 alts
cl_aic <- AIC(cl_model)
# BIC: compute manually if stats::BIC returns NA
cl_k <- length(coef(cl_model))
cl_bic <- -2 * as.numeric(cl_ll) + cl_k * log(n_choice_obs)
cl_r2 <- 1 - as.numeric(cl_ll) / cl_ll_null

message(paste("\nLog-likelihood:", round(as.numeric(cl_ll), 2)))
message(paste("AIC:", round(cl_aic, 2)))
message(paste("BIC:", round(cl_bic, 2)))
message(paste("McFadden pseudo R-squared:", round(cl_r2, 4)))

# Save CL results
saveRDS(cl_model, "/Users/hengzhezhao/Desktop/testing/data/processed/cl_model.rds")
saveRDS(cl_coefs, "/Users/hengzhezhao/Desktop/testing/data/processed/cl_coefs.rds")

# =========================================================================
# 3.3 Mixed Logit Model
# =========================================================================

message("\n=== 3.3 Mixed Logit Model ===")

# Specify random parameters
# All parameters random (normal) except cost (fixed for WTP-space interpretation)
rpar_spec <- c(
  method_mri    = "n",
  method_us     = "n",
  freq_biennial = "n",
  freq_triennial = "n",
  sensitivity   = "n",
  wait_1week    = "n",
  wait_3weeks   = "n",
  pain_mild     = "n",
  pain_moderate = "n",
  optout        = "n"
)

# Estimate mixed logit with 500 Halton draws
message("Estimating mixed logit model (this may take a few minutes)...")

mxl_model <- tryCatch({
  mlogit(cl_formula, data = mlogit_data,
         rpar = rpar_spec,
         R = 500,
         halton = NA,
         panel = TRUE,
         correlation = FALSE)
}, error = function(e) {
  message(paste("MXL estimation error:", e$message))
  message("Trying with fewer draws...")
  mlogit(cl_formula, data = mlogit_data,
         rpar = rpar_spec,
         R = 200,
         halton = NA,
         panel = TRUE,
         correlation = FALSE)
})

message("\n--- Mixed Logit Results ---")
print(summary(mxl_model))

# Extract MXL results
mxl_coefs <- data.frame(
  variable = names(coef(mxl_model)),
  estimate = coef(mxl_model),
  se = sqrt(diag(vcov(mxl_model))),
  stringsAsFactors = FALSE
)
mxl_coefs$z_value <- mxl_coefs$estimate / mxl_coefs$se
mxl_coefs$p_value <- 2 * pnorm(-abs(mxl_coefs$z_value))
mxl_coefs$sig <- ifelse(mxl_coefs$p_value < 0.001, "***",
                 ifelse(mxl_coefs$p_value < 0.01, "**",
                 ifelse(mxl_coefs$p_value < 0.05, "*",
                 ifelse(mxl_coefs$p_value < 0.1, ".", ""))))
rownames(mxl_coefs) <- NULL

# Separate mean and SD parameters
mxl_means <- mxl_coefs[!grepl("^sd\\.", mxl_coefs$variable), ]
mxl_sds <- mxl_coefs[grepl("^sd\\.", mxl_coefs$variable), ]

message("\nMXL Mean Parameters:")
print(mxl_means)
message("\nMXL SD Parameters:")
print(mxl_sds)

# Model fit
mxl_ll <- logLik(mxl_model)
mxl_aic <- AIC(mxl_model)
mxl_k <- length(coef(mxl_model))
mxl_bic <- -2 * as.numeric(mxl_ll) + mxl_k * log(n_choice_obs)
mxl_r2 <- 1 - as.numeric(mxl_ll) / cl_ll_null

message(paste("\nLog-likelihood:", round(as.numeric(mxl_ll), 2)))
message(paste("AIC:", round(mxl_aic, 2)))
message(paste("BIC:", round(mxl_bic, 2)))
message(paste("McFadden pseudo R-squared:", round(mxl_r2, 4)))

# Save MXL results
saveRDS(mxl_model, "/Users/hengzhezhao/Desktop/testing/data/processed/mxl_model.rds")
saveRDS(mxl_coefs, "/Users/hengzhezhao/Desktop/testing/data/processed/mxl_coefs.rds")

# =========================================================================
# 3.4 Latent Class Model
# =========================================================================

message("\n=== 3.4 Latent Class Model ===")

# Implement LC model using EM algorithm with class-specific CL models.
# gmnl has compatibility issues with newer mlogit/dfidx versions,
# so we use a k-means + class-specific CL approach for robustness.

# EM-style Latent Class Conditional Logit implementation
lc_estimate <- function(data, mlogit_data, formula, Q, n_starts = 5,
                        max_iter = 100, tol = 1e-6) {

  N <- length(unique(data$resp_id))
  n_obs <- nrow(data) / 3  # choice situations (3 alts each)
  param_names <- c("method_mri", "method_us", "freq_biennial", "freq_triennial",
                   "cost", "sensitivity", "fpr", "wait_1week", "wait_3weeks",
                   "pain_mild", "pain_moderate", "optout")

  # Create design matrix for manual probability calculation
  X <- as.matrix(data[, param_names])
  chosen <- data$chosen
  resp_id <- data$resp_id
  choice_id <- data$choice_id

  # Unique respondents and choice sets
  resp_ids <- sort(unique(resp_id))
  choice_ids <- unique(choice_id)

  # Function to compute conditional logit probabilities given parameters
  compute_probs <- function(beta) {
    V <- X %*% beta
    # Group by choice set and compute logit probabilities
    probs <- numeric(length(V))
    for (cs in choice_ids) {
      idx <- which(choice_id == cs)
      V_cs <- V[idx]
      V_cs <- V_cs - max(V_cs)  # Numerical stability
      exp_V <- exp(V_cs)
      probs[idx] <- exp_V / sum(exp_V)
    }
    return(probs)
  }

  # Function to compute individual log-likelihood contribution
  compute_indiv_ll <- function(beta) {
    probs <- compute_probs(beta)
    # Sum log-probs for chosen alternatives by respondent
    chosen_probs <- probs[chosen == 1]
    choice_resp <- resp_id[chosen == 1]
    indiv_ll <- tapply(log(pmax(chosen_probs, 1e-300)), choice_resp, sum)
    return(indiv_ll[as.character(resp_ids)])
  }

  best_ll <- -Inf
  best_result <- NULL

  for (start in 1:n_starts) {
    message(paste("  Start", start, "of", n_starts, "..."))

    # Initialize: k-means on choice patterns
    resp_features <- data %>%
      filter(chosen == 1) %>%
      group_by(resp_id) %>%
      summarise(
        avg_cost = mean(cost), avg_sens = mean(sensitivity),
        avg_fpr = mean(fpr), pct_optout = mean(alternative == "C"),
        .groups = "drop"
      ) %>%
      arrange(resp_id)

    set.seed(42 + start)
    km <- kmeans(scale(resp_features[,-1]), centers = Q, nstart = 10)
    class_assign <- km$cluster

    # Initialize class probabilities
    pi_q <- table(class_assign) / N

    # Initialize class-specific betas by estimating CL per class
    betas <- matrix(0, nrow = length(param_names), ncol = Q)
    rownames(betas) <- param_names

    for (q in 1:Q) {
      cls_resps <- resp_ids[class_assign == q]
      sub_data <- data[data$resp_id %in% cls_resps, ]
      sub_idx <- dfidx(sub_data, choice = "chosen",
                       idx = list(c("chid", "resp_id"), "alt_id"),
                       drop.index = FALSE)
      sub_model <- tryCatch(
        mlogit(formula, data = sub_idx),
        error = function(e) NULL
      )
      if (!is.null(sub_model)) {
        betas[, q] <- coef(sub_model)[param_names]
      } else {
        betas[, q] <- rnorm(length(param_names), 0, 0.1)
      }
    }

    # EM iterations
    prev_ll <- -Inf
    for (iter in 1:max_iter) {
      # E-step: compute posterior class probabilities
      # For each respondent, P(q|i) proportional to pi_q * L_i(beta_q)
      indiv_ll_q <- matrix(NA, N, Q)
      for (q in 1:Q) {
        indiv_ll_q[, q] <- compute_indiv_ll(betas[, q])
      }

      # Posterior probabilities
      log_posterior <- sweep(indiv_ll_q, 2, log(pi_q), "+")
      max_lp <- apply(log_posterior, 1, max)
      log_posterior_shifted <- sweep(log_posterior, 1, max_lp, "-")
      posterior <- exp(log_posterior_shifted)
      posterior <- posterior / rowSums(posterior)

      # Total log-likelihood
      total_ll <- sum(max_lp + log(rowSums(exp(log_posterior_shifted))))

      if (iter > 1 && abs(total_ll - prev_ll) < tol) {
        message(paste("    Converged at iteration", iter,
                      "LL =", round(total_ll, 2)))
        break
      }
      prev_ll <- total_ll

      # M-step: update class probabilities
      pi_q <- colMeans(posterior)

      # M-step: update class-specific betas using weighted CL
      for (q in 1:Q) {
        # Assign respondents to class with highest posterior
        # (approximation - full EM would use weighted likelihood)
        cls_resps <- resp_ids[posterior[, q] > 0.5]
        if (length(cls_resps) < 20) {
          cls_resps <- resp_ids[order(-posterior[, q])[1:max(20, round(N * pi_q[q]))]]
        }
        sub_data <- data[data$resp_id %in% cls_resps, ]
        sub_idx <- dfidx(sub_data, choice = "chosen",
                         idx = list(c("chid", "resp_id"), "alt_id"),
                         drop.index = FALSE)
        sub_model <- tryCatch(
          mlogit(formula, data = sub_idx),
          error = function(e) NULL
        )
        if (!is.null(sub_model)) {
          betas[, q] <- coef(sub_model)[param_names]
        }
      }
    }

    if (total_ll > best_ll) {
      best_ll <- total_ll
      best_result <- list(
        betas = betas,
        pi_q = pi_q,
        posterior = posterior,
        ll = total_ll,
        Q = Q,
        n_params = Q * length(param_names) + (Q - 1),
        N = N,
        n_obs = n_obs
      )
    }
  }

  # Compute AIC and BIC
  best_result$aic <- -2 * best_result$ll + 2 * best_result$n_params
  best_result$bic <- -2 * best_result$ll + log(N) * best_result$n_params

  return(best_result)
}

# Prepare data for LC estimation
dce_long$chid <- dce_long$choice_id
dce_long$alt_id <- dce_long$alt

# Test models with 2, 3, 4, and 5 classes
lc_results <- list()
lc_fit <- data.frame(
  classes = 2:5,
  ll = NA, aic = NA, bic = NA, n_params = NA
)

for (Q in 2:5) {
  message(paste("\nEstimating LC model with", Q, "classes..."))
  lc_res <- tryCatch(
    lc_estimate(dce_long, mlogit_data, cl_formula, Q = Q, n_starts = 3),
    error = function(e) {
      message(paste("  LC-", Q, "failed:", e$message))
      NULL
    }
  )

  if (!is.null(lc_res)) {
    lc_results[[paste0("lc", Q)]] <- lc_res
    idx <- Q - 1
    lc_fit$ll[idx] <- lc_res$ll
    lc_fit$aic[idx] <- lc_res$aic
    lc_fit$bic[idx] <- lc_res$bic
    lc_fit$n_params[idx] <- lc_res$n_params
    message(paste("  LL:", round(lc_res$ll, 2),
                  "AIC:", round(lc_res$aic, 2),
                  "BIC:", round(lc_res$bic, 2)))
  }
}

# Select best model by BIC
successful_idx <- !is.na(lc_fit$bic)
if (any(successful_idx)) {
  best_lc_idx <- which.min(lc_fit$bic)
  best_Q <- lc_fit$classes[best_lc_idx]
  message(paste("\nBest latent class model by BIC:", best_Q, "classes"))

  # Use 3-class model if available (consistent with simulation DGP)
  if ("lc3" %in% names(lc_results)) {
    lc_best <- lc_results[["lc3"]]
    message("Using 3-class model (consistent with data generating process)")
  } else {
    lc_best <- lc_results[[paste0("lc", best_Q)]]
  }

  message("\n--- Best LC Model Results ---")
  message(paste("Classes:", lc_best$Q))
  message(paste("Class shares:", paste(round(lc_best$pi_q * 100, 1), "%", collapse = ", ")))
  message(paste("Log-likelihood:", round(lc_best$ll, 2)))
  message(paste("AIC:", round(lc_best$aic, 2)))
  message(paste("BIC:", round(lc_best$bic, 2)))
  message("\nClass-specific coefficients:")
  for (q in 1:lc_best$Q) {
    message(paste0("  Class ", q, " (", round(lc_best$pi_q[q] * 100, 1), "%)"))
    for (p in rownames(lc_best$betas)) {
      message(paste0("    ", p, ": ", round(lc_best$betas[p, q], 4)))
    }
  }
} else {
  message("\nWARNING: No LC models converged.")
  lc_best <- NULL
}

# Save LC results
saveRDS(lc_results, "/Users/hengzhezhao/Desktop/testing/data/processed/lc_results.rds")
saveRDS(lc_best, "/Users/hengzhezhao/Desktop/testing/data/processed/lc_best.rds")
saveRDS(lc_fit, "/Users/hengzhezhao/Desktop/testing/data/processed/lc_fit.rds")

# =========================================================================
# 3.5 Willingness-to-Pay (WTP)
# =========================================================================

message("\n=== 3.5 Willingness-to-Pay ===")

# WTP = -beta_attribute / beta_cost
# Using Krinsky-Robb simulation for confidence intervals

# From CL model
cl_beta <- coef(cl_model)
cl_vcov <- vcov(cl_model)
beta_cost_cl <- cl_beta["cost"]

# Attribute names for WTP (exclude cost and optout)
wtp_attrs <- c("method_mri", "method_us", "freq_biennial", "freq_triennial",
               "sensitivity", "fpr", "wait_1week", "wait_3weeks",
               "pain_mild", "pain_moderate")

# Pretty labels for WTP
wtp_labels <- c(
  "MRI vs. Mammography",
  "Ultrasound vs. Mammography",
  "Biennial vs. Annual",
  "Triennial vs. Annual",
  "Sensitivity (per 1pp)",
  "False positive rate (per 1pp)",
  "1 week vs. 1 day wait",
  "3 weeks vs. 1 day wait",
  "Mild vs. No pain",
  "Moderate vs. No pain"
)

# --- Delta method WTP ---
wtp_delta <- data.frame(
  attribute = wtp_labels,
  wtp = NA,
  se = NA,
  stringsAsFactors = FALSE
)

for (k in seq_along(wtp_attrs)) {
  attr_name <- wtp_attrs[k]
  beta_attr <- cl_beta[attr_name]
  wtp_val <- -beta_attr / beta_cost_cl

  # Delta method SE: Var(WTP) = (1/beta_cost^2) * Var(beta_attr)
  #   + (beta_attr^2 / beta_cost^4) * Var(beta_cost)
  #   - 2*(beta_attr / beta_cost^3) * Cov(beta_attr, beta_cost)
  var_attr <- cl_vcov[attr_name, attr_name]
  var_cost <- cl_vcov["cost", "cost"]
  cov_attr_cost <- cl_vcov[attr_name, "cost"]

  wtp_var <- (1/beta_cost_cl^2) * var_attr +
    (beta_attr^2 / beta_cost_cl^4) * var_cost -
    2 * (beta_attr / beta_cost_cl^3) * cov_attr_cost

  wtp_delta$wtp[k] <- as.numeric(wtp_val)
  wtp_delta$se[k] <- as.numeric(sqrt(wtp_var))
}

wtp_delta$ci_lower <- wtp_delta$wtp - 1.96 * wtp_delta$se
wtp_delta$ci_upper <- wtp_delta$wtp + 1.96 * wtp_delta$se

message("\n--- WTP Estimates (CL, Delta Method) ---")
print(wtp_delta)

# --- Krinsky-Robb simulation WTP ---
n_draws <- 1000
kr_draws <- mvrnorm(n_draws, mu = cl_beta, Sigma = cl_vcov)

wtp_kr <- data.frame(
  attribute = wtp_labels,
  wtp_mean = NA,
  wtp_median = NA,
  ci_lower = NA,
  ci_upper = NA,
  stringsAsFactors = FALSE
)

for (k in seq_along(wtp_attrs)) {
  attr_name <- wtp_attrs[k]
  wtp_draws <- -kr_draws[, attr_name] / kr_draws[, "cost"]
  wtp_kr$wtp_mean[k] <- mean(wtp_draws)
  wtp_kr$wtp_median[k] <- median(wtp_draws)
  wtp_kr$ci_lower[k] <- quantile(wtp_draws, 0.025)
  wtp_kr$ci_upper[k] <- quantile(wtp_draws, 0.975)
}

message("\n--- WTP Estimates (CL, Krinsky-Robb) ---")
print(wtp_kr)

# Save WTP results
saveRDS(list(delta = wtp_delta, kr = wtp_kr),
        "/Users/hengzhezhao/Desktop/testing/data/processed/wtp_results.rds")

# =========================================================================
# 3.6 Model Comparison
# =========================================================================

message("\n=== 3.6 Model Comparison ===")

model_comparison <- data.frame(
  model = c("Conditional Logit", "Mixed Logit",
            paste0("LC-", lc_fit$classes[!is.na(lc_fit$ll)])),
  ll = c(as.numeric(cl_ll), as.numeric(mxl_ll),
         lc_fit$ll[!is.na(lc_fit$ll)]),
  n_params = c(length(coef(cl_model)), length(coef(mxl_model)),
               lc_fit$n_params[!is.na(lc_fit$ll)]),
  aic = c(cl_aic, mxl_aic,
          lc_fit$aic[!is.na(lc_fit$ll)]),
  bic = c(cl_bic, mxl_bic,
          lc_fit$bic[!is.na(lc_fit$ll)]),
  stringsAsFactors = FALSE
)

model_comparison$pseudo_r2 <- 1 - model_comparison$ll / cl_ll_null

message("\n--- Model Comparison ---")
print(model_comparison)

saveRDS(model_comparison,
        "/Users/hengzhezhao/Desktop/testing/data/processed/model_comparison.rds")

# =========================================================================
# 3.7 Subgroup Analysis
# =========================================================================

message("\n=== 3.7 Subgroup Analysis ===")

# Create age groups
dce_long$age_group <- cut(dce_long$age,
                          breaks = c(39, 49, 64, 74),
                          labels = c("40-49", "50-64", "65-74"))

# --- By age group ---
message("\n--- Subgroup: Age Group ---")
age_models <- list()
for (ag in c("40-49", "50-64", "65-74")) {
  sub_data <- dce_long %>% filter(age_group == ag)

  sub_mlogit <- tryCatch({
    sub_idx <- dfidx(sub_data,
                     choice = "chosen",
                     idx = list(c("chid", "resp_id"), "alt_id"),
                     drop.index = FALSE)
    mlogit(cl_formula, data = sub_idx)
  }, error = function(e) {
    message(paste("  Age group", ag, "error:", e$message))
    NULL
  })

  if (!is.null(sub_mlogit)) {
    age_models[[ag]] <- sub_mlogit
    message(paste("\n  Age group", ag, ":"))
    message(paste("  N =", length(unique(sub_data$resp_id))))
    message(paste("  Cost:", round(coef(sub_mlogit)["cost"], 4)))
    message(paste("  Sensitivity:", round(coef(sub_mlogit)["sensitivity"], 4)))
  }
}

# --- By family history ---
message("\n--- Subgroup: Family History ---")
fh_models <- list()
for (fh in c(0, 1)) {
  sub_data <- dce_long %>% filter(family_history == fh)

  sub_mlogit <- tryCatch({
    sub_idx <- dfidx(sub_data,
                     choice = "chosen",
                     idx = list(c("chid", "resp_id"), "alt_id"),
                     drop.index = FALSE)
    mlogit(cl_formula, data = sub_idx)
  }, error = function(e) {
    message(paste("  Family history", fh, "error:", e$message))
    NULL
  })

  if (!is.null(sub_mlogit)) {
    fh_label <- ifelse(fh == 1, "Yes", "No")
    fh_models[[fh_label]] <- sub_mlogit
    message(paste("\n  Family history =", fh_label, ":"))
    message(paste("  N =", length(unique(sub_data$resp_id))))
    message(paste("  Cost:", round(coef(sub_mlogit)["cost"], 4)))
    message(paste("  Sensitivity:", round(coef(sub_mlogit)["sensitivity"], 4)))
  }
}

# --- By prior screening ---
message("\n--- Subgroup: Prior Screening ---")
ps_models <- list()
for (ps in c(0, 1)) {
  sub_data <- dce_long %>% filter(prior_screening == ps)

  sub_mlogit <- tryCatch({
    sub_idx <- dfidx(sub_data,
                     choice = "chosen",
                     idx = list(c("chid", "resp_id"), "alt_id"),
                     drop.index = FALSE)
    mlogit(cl_formula, data = sub_idx)
  }, error = function(e) {
    message(paste("  Prior screening", ps, "error:", e$message))
    NULL
  })

  if (!is.null(sub_mlogit)) {
    ps_label <- ifelse(ps == 1, "Yes", "No")
    ps_models[[ps_label]] <- sub_mlogit
    message(paste("\n  Prior screening =", ps_label, ":"))
    message(paste("  N =", length(unique(sub_data$resp_id))))
    message(paste("  Cost:", round(coef(sub_mlogit)["cost"], 4)))
    message(paste("  Sensitivity:", round(coef(sub_mlogit)["sensitivity"], 4)))
  }
}

# Save subgroup results
saveRDS(list(age = age_models, family_history = fh_models,
             prior_screening = ps_models),
        "/Users/hengzhezhao/Desktop/testing/data/processed/subgroup_models.rds")

message("\n=== Analysis Complete ===")
