###############################################################################
# 05_tables.R
# Generate LaTeX tables for DCE analysis results
# Project: Health Preferences for Breast Cancer Screening
###############################################################################

library(mlogit)
library(dplyr)
library(tidyr)
library(xtable)
library(kableExtra)

set.seed(42)

message("=== Generating LaTeX Tables ===")

# -------------------------------------------------------------------------
# Load results
# -------------------------------------------------------------------------

respondents <- read.csv("/Users/hengzhezhao/Desktop/testing/data/processed/respondents.csv",
                        stringsAsFactors = TRUE)
cl_coefs <- readRDS("/Users/hengzhezhao/Desktop/testing/data/processed/cl_coefs.rds")
cl_model <- readRDS("/Users/hengzhezhao/Desktop/testing/data/processed/cl_model.rds")
mxl_coefs <- readRDS("/Users/hengzhezhao/Desktop/testing/data/processed/mxl_coefs.rds")
mxl_model <- readRDS("/Users/hengzhezhao/Desktop/testing/data/processed/mxl_model.rds")
lc_best <- readRDS("/Users/hengzhezhao/Desktop/testing/data/processed/lc_best.rds")
lc_fit <- readRDS("/Users/hengzhezhao/Desktop/testing/data/processed/lc_fit.rds")
wtp_results <- readRDS("/Users/hengzhezhao/Desktop/testing/data/processed/wtp_results.rds")
model_comparison <- readRDS("/Users/hengzhezhao/Desktop/testing/data/processed/model_comparison.rds")
demo_tables <- readRDS("/Users/hengzhezhao/Desktop/testing/data/processed/demo_tables.rds")
choice_summary <- readRDS("/Users/hengzhezhao/Desktop/testing/data/processed/choice_summary.rds")

# =========================================================================
# Table 1: Sample Demographics
# =========================================================================

message("Generating Table 1: Demographics...")

N <- nrow(respondents)

# Build demographics table
demo_rows <- data.frame(
  Variable = character(),
  Category = character(),
  Value = character(),
  stringsAsFactors = FALSE
)

# Age
demo_rows <- rbind(demo_rows, data.frame(
  Variable = "Age (years)",
  Category = "Mean (SD)",
  Value = paste0(round(mean(respondents$age), 1), " (", round(sd(respondents$age), 1), ")"),
  stringsAsFactors = FALSE
))
demo_rows <- rbind(demo_rows, data.frame(
  Variable = "",
  Category = "Range",
  Value = paste0(min(respondents$age), "--", max(respondents$age)),
  stringsAsFactors = FALSE
))

# Income
for (lvl in levels(respondents$income)) {
  n_lvl <- sum(respondents$income == lvl)
  pct <- round(n_lvl / N * 100, 1)
  demo_rows <- rbind(demo_rows, data.frame(
    Variable = ifelse(lvl == levels(respondents$income)[1], "Income", ""),
    Category = lvl,
    Value = paste0(n_lvl, " (", pct, "\\%)"),
    stringsAsFactors = FALSE
  ))
}

# Education
for (lvl in levels(respondents$education)) {
  n_lvl <- sum(respondents$education == lvl)
  pct <- round(n_lvl / N * 100, 1)
  demo_rows <- rbind(demo_rows, data.frame(
    Variable = ifelse(lvl == levels(respondents$education)[1], "Education", ""),
    Category = lvl,
    Value = paste0(n_lvl, " (", pct, "\\%)"),
    stringsAsFactors = FALSE
  ))
}

# Family history
n_fh <- sum(respondents$family_history == 1)
demo_rows <- rbind(demo_rows, data.frame(
  Variable = "Family history",
  Category = "Yes",
  Value = paste0(n_fh, " (", round(n_fh/N*100, 1), "\\%)"),
  stringsAsFactors = FALSE
))

# Prior screening
n_ps <- sum(respondents$prior_screening == 1)
demo_rows <- rbind(demo_rows, data.frame(
  Variable = "Prior screening",
  Category = "Yes",
  Value = paste0(n_ps, " (", round(n_ps/N*100, 1), "\\%)"),
  stringsAsFactors = FALSE
))

# Insurance
for (lvl in levels(respondents$insurance)) {
  n_lvl <- sum(respondents$insurance == lvl)
  pct <- round(n_lvl / N * 100, 1)
  demo_rows <- rbind(demo_rows, data.frame(
    Variable = ifelse(lvl == levels(respondents$insurance)[1], "Insurance", ""),
    Category = lvl,
    Value = paste0(n_lvl, " (", pct, "\\%)"),
    stringsAsFactors = FALSE
  ))
}

# Race
for (lvl in levels(respondents$race)) {
  n_lvl <- sum(respondents$race == lvl)
  pct <- round(n_lvl / N * 100, 1)
  demo_rows <- rbind(demo_rows, data.frame(
    Variable = ifelse(lvl == levels(respondents$race)[1], "Race/Ethnicity", ""),
    Category = lvl,
    Value = paste0(n_lvl, " (", pct, "\\%)"),
    stringsAsFactors = FALSE
  ))
}

# Write LaTeX table
tab1_tex <- paste0(
  "\\begin{table}[htbp]\n",
  "\\centering\n",
  "\\caption{Sample Characteristics (N = ", N, ")}\n",
  "\\label{tab:demographics}\n",
  "\\begin{tabular}{llr}\n",
  "\\toprule\n",
  "Variable & Category & n (\\%) or Mean (SD) \\\\\n",
  "\\midrule\n"
)

for (i in 1:nrow(demo_rows)) {
  tab1_tex <- paste0(tab1_tex,
                     demo_rows$Variable[i], " & ",
                     demo_rows$Category[i], " & ",
                     demo_rows$Value[i], " \\\\\n")
  # Add midrule between variable groups
  if (i < nrow(demo_rows) && demo_rows$Variable[i+1] != "" && i > 1) {
    tab1_tex <- paste0(tab1_tex, "\\midrule\n")
  }
}

tab1_tex <- paste0(tab1_tex,
                   "\\bottomrule\n",
                   "\\end{tabular}\n",
                   "\\end{table}\n")

writeLines(tab1_tex, "/Users/hengzhezhao/Desktop/testing/tables/tab_demographics.tex")
message("  Saved tab_demographics.tex")

# =========================================================================
# Table 2: DCE Attributes and Levels
# =========================================================================

message("Generating Table 2: Attributes and levels...")

tab2_tex <- paste0(
  "\\begin{table}[htbp]\n",
  "\\centering\n",
  "\\caption{DCE Attributes and Levels}\n",
  "\\label{tab:attributes}\n",
  "\\begin{tabular}{lll}\n",
  "\\toprule\n",
  "Attribute & Levels & Coding \\\\\n",
  "\\midrule\n",
  "Screening method & Mammography (ref.), MRI, Ultrasound & Dummy \\\\\n",
  "Screening frequency & Annual (ref.), Biennial, Triennial & Dummy \\\\\n",
  "Out-of-pocket cost & \\$0, \\$50, \\$150, \\$300 & Continuous \\\\\n",
  "Test sensitivity & 70\\%, 85\\%, 95\\% & Continuous \\\\\n",
  "False positive rate & 5\\%, 10\\%, 15\\% & Continuous \\\\\n",
  "Waiting time for results & 1 day (ref.), 1 week, 3 weeks & Dummy \\\\\n",
  "Pain/discomfort & None (ref.), Mild, Moderate & Dummy \\\\\n",
  "\\bottomrule\n",
  "\\end{tabular}\n",
  "\\begin{tablenotes}\n",
  "\\small\n",
  "\\item Note: ref.\\ = reference level. Continuous variables entered as linear terms.\n",
  "\\item Design: D-optimal fractional factorial, 36 choice sets in 3 blocks of 12.\n",
  "\\item Each choice set: 2 screening alternatives + opt-out (no screening).\n",
  "\\end{tablenotes}\n",
  "\\end{table}\n"
)

writeLines(tab2_tex, "/Users/hengzhezhao/Desktop/testing/tables/tab_attributes.tex")
message("  Saved tab_attributes.tex")

# =========================================================================
# Table 3: Conditional Logit Results
# =========================================================================

message("Generating Table 3: Conditional logit results...")

# Pretty variable labels
var_labels <- c(
  "method_mri" = "MRI vs.\\ Mammography",
  "method_us" = "Ultrasound vs.\\ Mammography",
  "freq_biennial" = "Biennial vs.\\ Annual",
  "freq_triennial" = "Triennial vs.\\ Annual",
  "cost" = "Cost (per \\$1)",
  "sensitivity" = "Sensitivity (per 1 pp)",
  "fpr" = "False positive rate (per 1 pp)",
  "wait_1week" = "1 week vs.\\ 1 day wait",
  "wait_3weeks" = "3 weeks vs.\\ 1 day wait",
  "pain_mild" = "Mild vs.\\ None",
  "pain_moderate" = "Moderate vs.\\ None",
  "optout" = "Opt-out (no screening)"
)

cl_tab <- cl_coefs %>%
  mutate(
    label = var_labels[variable],
    coef_str = paste0(sprintf("%.4f", estimate), sig),
    se_str = paste0("(", sprintf("%.4f", se), ")"),
    or = exp(estimate),
    or_str = sprintf("%.3f", or)
  )

# Build LaTeX
tab3_tex <- paste0(
  "\\begin{table}[htbp]\n",
  "\\centering\n",
  "\\caption{Conditional Logit Model Results}\n",
  "\\label{tab:cl_results}\n",
  "\\begin{tabular}{lrrr}\n",
  "\\toprule\n",
  "Attribute Level & Coefficient (SE) & $p$-value & OR \\\\\n",
  "\\midrule\n"
)

for (i in 1:nrow(cl_tab)) {
  tab3_tex <- paste0(tab3_tex,
                     cl_tab$label[i], " & ",
                     cl_tab$coef_str[i], " & ",
                     sprintf("%.4f", cl_tab$p_value[i]), " & ",
                     cl_tab$or_str[i], " \\\\\n",
                     " & ", cl_tab$se_str[i], " & & \\\\\n")
}

# Model fit
cl_ll <- as.numeric(logLik(cl_model))
n_long <- nrow(read.csv("/Users/hengzhezhao/Desktop/testing/data/processed/dce_long.csv"))
n_choice_obs <- n_long / 3
cl_ll_null <- -log(3) * n_choice_obs
cl_k <- length(coef(cl_model))
cl_bic_val <- -2 * cl_ll + cl_k * log(n_choice_obs)
cl_r2 <- 1 - cl_ll / cl_ll_null

tab3_tex <- paste0(tab3_tex,
                   "\\midrule\n",
                   "Log-likelihood & \\multicolumn{3}{c}{", sprintf("%.2f", cl_ll), "} \\\\\n",
                   "AIC & \\multicolumn{3}{c}{", sprintf("%.2f", AIC(cl_model)), "} \\\\\n",
                   "BIC & \\multicolumn{3}{c}{", sprintf("%.2f", cl_bic_val), "} \\\\\n",
                   "McFadden $R^2$ & \\multicolumn{3}{c}{", sprintf("%.4f", cl_r2), "} \\\\\n",
                   "N (choice obs.) & \\multicolumn{3}{c}{", n_choice_obs, "} \\\\\n",
                   "\\bottomrule\n",
                   "\\end{tabular}\n",
                   "\\begin{tablenotes}\n",
                   "\\small\n",
                   "\\item Note: *** $p < 0.001$, ** $p < 0.01$, * $p < 0.05$, . $p < 0.1$.\n",
                   "\\item OR = odds ratio. SE in parentheses.\n",
                   "\\end{tablenotes}\n",
                   "\\end{table}\n")

writeLines(tab3_tex, "/Users/hengzhezhao/Desktop/testing/tables/tab_cl_results.tex")
message("  Saved tab_cl_results.tex")

# =========================================================================
# Table 4: Mixed Logit Results
# =========================================================================

message("Generating Table 4: Mixed logit results...")

mxl_means <- mxl_coefs %>% filter(!grepl("^sd\\.", variable))
mxl_sds <- mxl_coefs %>% filter(grepl("^sd\\.", variable))

# Match means and SDs
mxl_tab <- mxl_means %>%
  mutate(label = ifelse(variable %in% names(var_labels),
                        var_labels[variable], variable))

# Find corresponding SDs
mxl_tab$sd_est <- NA
mxl_tab$sd_se <- NA
mxl_tab$sd_sig <- NA
for (i in 1:nrow(mxl_tab)) {
  sd_name <- paste0("sd.", mxl_tab$variable[i])
  sd_row <- mxl_sds[mxl_sds$variable == sd_name, ]
  if (nrow(sd_row) == 1) {
    mxl_tab$sd_est[i] <- sd_row$estimate
    mxl_tab$sd_se[i] <- sd_row$se
    mxl_tab$sd_sig[i] <- sd_row$sig
  }
}

tab4_tex <- paste0(
  "\\begin{table}[htbp]\n",
  "\\centering\n",
  "\\caption{Mixed Logit Model Results}\n",
  "\\label{tab:mxl_results}\n",
  "\\begin{tabular}{lrrrr}\n",
  "\\toprule\n",
  " & \\multicolumn{2}{c}{Mean} & \\multicolumn{2}{c}{SD} \\\\\n",
  "\\cmidrule(lr){2-3} \\cmidrule(lr){4-5}\n",
  "Attribute Level & Coef. & SE & Coef. & SE \\\\\n",
  "\\midrule\n"
)

for (i in 1:nrow(mxl_tab)) {
  sd_str <- ifelse(is.na(mxl_tab$sd_est[i]), "--",
                   paste0(sprintf("%.4f", mxl_tab$sd_est[i]), mxl_tab$sd_sig[i]))
  sd_se_str <- ifelse(is.na(mxl_tab$sd_se[i]), "--",
                      sprintf("%.4f", mxl_tab$sd_se[i]))

  tab4_tex <- paste0(tab4_tex,
                     mxl_tab$label[i], " & ",
                     sprintf("%.4f", mxl_tab$estimate[i]), mxl_tab$sig[i], " & ",
                     sprintf("%.4f", mxl_tab$se[i]), " & ",
                     sd_str, " & ",
                     sd_se_str, " \\\\\n")
}

mxl_ll <- as.numeric(logLik(mxl_model))
mxl_k <- length(coef(mxl_model))
mxl_bic_val <- -2 * mxl_ll + mxl_k * log(n_choice_obs)
mxl_r2 <- 1 - mxl_ll / cl_ll_null

tab4_tex <- paste0(tab4_tex,
                   "\\midrule\n",
                   "Log-likelihood & \\multicolumn{4}{c}{", sprintf("%.2f", mxl_ll), "} \\\\\n",
                   "AIC & \\multicolumn{4}{c}{", sprintf("%.2f", AIC(mxl_model)), "} \\\\\n",
                   "BIC & \\multicolumn{4}{c}{", sprintf("%.2f", mxl_bic_val), "} \\\\\n",
                   "McFadden $R^2$ & \\multicolumn{4}{c}{", sprintf("%.4f", mxl_r2), "} \\\\\n",
                   "Halton draws & \\multicolumn{4}{c}{500} \\\\\n",
                   "\\bottomrule\n",
                   "\\end{tabular}\n",
                   "\\begin{tablenotes}\n",
                   "\\small\n",
                   "\\item Note: *** $p < 0.001$, ** $p < 0.01$, * $p < 0.05$, . $p < 0.1$.\n",
                   "\\item All parameters except cost and FPR specified as random (normal distribution).\n",
                   "\\item SD = standard deviation of random parameter distribution.\n",
                   "\\end{tablenotes}\n",
                   "\\end{table}\n")

writeLines(tab4_tex, "/Users/hengzhezhao/Desktop/testing/tables/tab_mxl_results.tex")
message("  Saved tab_mxl_results.tex")

# =========================================================================
# Table 5: Latent Class Results
# =========================================================================

message("Generating Table 5: Latent class results...")

# LC results use custom format: lc_best$betas, lc_best$pi_q, lc_best$Q, lc_best$ll, etc.
if (!is.null(lc_best) && !is.null(lc_best$betas)) {
  Q <- lc_best$Q
  param_names <- rownames(lc_best$betas)

  tab5_tex <- paste0(
    "\\begin{table}[htbp]\n",
    "\\centering\n",
    "\\caption{Latent Class Model Results (", Q, " Classes)}\n",
    "\\label{tab:lc_results}\n",
    "\\begin{tabular}{l", paste(rep("r", Q), collapse = ""), "}\n",
    "\\toprule\n",
    " & ", paste(paste0("Class ", 1:Q, " (", round(lc_best$pi_q * 100, 1), "\\%)"),
                 collapse = " & "), " \\\\\n",
    "\\midrule\n"
  )

  for (pname in param_names) {
    label <- ifelse(pname %in% names(var_labels), var_labels[pname], pname)
    vals <- character(Q)
    for (q in 1:Q) {
      vals[q] <- sprintf("%.4f", lc_best$betas[pname, q])
    }
    tab5_tex <- paste0(tab5_tex, label, " & ",
                       paste(vals, collapse = " & "), " \\\\\n")
  }

  n_choice_obs <- nrow(read.csv("/Users/hengzhezhao/Desktop/testing/data/processed/dce_long.csv")) / 3
  lc_ll_null <- -log(3) * n_choice_obs
  lc_ll <- lc_best$ll
  lc_r2 <- 1 - lc_ll / lc_ll_null

  tab5_tex <- paste0(tab5_tex,
                     "\\midrule\n",
                     "Class share (\\%) & ",
                     paste(sprintf("%.1f", lc_best$pi_q * 100), collapse = " & "),
                     " \\\\\n",
                     "\\midrule\n",
                     "Log-likelihood & \\multicolumn{", Q, "}{c}{", sprintf("%.2f", lc_ll), "} \\\\\n",
                     "AIC & \\multicolumn{", Q, "}{c}{", sprintf("%.2f", lc_best$aic), "} \\\\\n",
                     "BIC & \\multicolumn{", Q, "}{c}{", sprintf("%.2f", lc_best$bic), "} \\\\\n",
                     "McFadden $R^2$ & \\multicolumn{", Q, "}{c}{", sprintf("%.4f", lc_r2), "} \\\\\n",
                     "\\bottomrule\n",
                     "\\end{tabular}\n",
                     "\\begin{tablenotes}\n",
                     "\\small\n",
                     "\\item Note: Class-specific coefficients estimated via EM algorithm with\n",
                     "\\item class-specific conditional logit models.\n",
                     "\\end{tablenotes}\n",
                     "\\end{table}\n")
} else {
  tab5_tex <- paste0(
    "\\begin{table}[htbp]\n",
    "\\centering\n",
    "\\caption{Latent Class Model Results}\n",
    "\\label{tab:lc_results}\n",
    "\\begin{tabular}{l}\n",
    "\\toprule\n",
    "Latent class model results not available. \\\\\n",
    "\\bottomrule\n",
    "\\end{tabular}\n",
    "\\end{table}\n"
  )
}

writeLines(tab5_tex, "/Users/hengzhezhao/Desktop/testing/tables/tab_lc_results.tex")
message("  Saved tab_lc_results.tex")

# =========================================================================
# Table 6: WTP Estimates
# =========================================================================

message("Generating Table 6: WTP estimates...")

wtp_kr <- wtp_results$kr

tab6_tex <- paste0(
  "\\begin{table}[htbp]\n",
  "\\centering\n",
  "\\caption{Willingness-to-Pay Estimates (Conditional Logit)}\n",
  "\\label{tab:wtp}\n",
  "\\begin{tabular}{lrrr}\n",
  "\\toprule\n",
  "Attribute Level Change & WTP (\\$) & 95\\% CI Lower & 95\\% CI Upper \\\\\n",
  "\\midrule\n"
)

for (i in 1:nrow(wtp_kr)) {
  tab6_tex <- paste0(tab6_tex,
                     wtp_kr$attribute[i], " & ",
                     sprintf("%.2f", wtp_kr$wtp_mean[i]), " & ",
                     sprintf("%.2f", wtp_kr$ci_lower[i]), " & ",
                     sprintf("%.2f", wtp_kr$ci_upper[i]), " \\\\\n")
}

tab6_tex <- paste0(tab6_tex,
                   "\\bottomrule\n",
                   "\\end{tabular}\n",
                   "\\begin{tablenotes}\n",
                   "\\small\n",
                   "\\item Note: WTP = $-\\hat{\\beta}_{\\text{attribute}} / \\hat{\\beta}_{\\text{cost}}$.\n",
                   "\\item 95\\% CI computed using Krinsky-Robb simulation (1,000 draws).\n",
                   "\\item Positive values indicate willingness to pay for the attribute level change.\n",
                   "\\item Negative values indicate compensation required to accept the attribute level.\n",
                   "\\end{tablenotes}\n",
                   "\\end{table}\n")

writeLines(tab6_tex, "/Users/hengzhezhao/Desktop/testing/tables/tab_wtp.tex")
message("  Saved tab_wtp.tex")

# =========================================================================
# Table 7: Model Comparison
# =========================================================================

message("Generating Table 7: Model comparison...")

tab7_tex <- paste0(
  "\\begin{table}[htbp]\n",
  "\\centering\n",
  "\\caption{Model Comparison}\n",
  "\\label{tab:model_comparison}\n",
  "\\begin{tabular}{lrrrrr}\n",
  "\\toprule\n",
  "Model & LL & $K$ & AIC & BIC & Pseudo-$R^2$ \\\\\n",
  "\\midrule\n"
)

for (i in 1:nrow(model_comparison)) {
  tab7_tex <- paste0(tab7_tex,
                     model_comparison$model[i], " & ",
                     sprintf("%.2f", model_comparison$ll[i]), " & ",
                     model_comparison$n_params[i], " & ",
                     sprintf("%.2f", model_comparison$aic[i]), " & ",
                     sprintf("%.2f", model_comparison$bic[i]), " & ",
                     sprintf("%.4f", model_comparison$pseudo_r2[i]), " \\\\\n")
}

tab7_tex <- paste0(tab7_tex,
                   "\\bottomrule\n",
                   "\\end{tabular}\n",
                   "\\begin{tablenotes}\n",
                   "\\small\n",
                   "\\item Note: LL = log-likelihood; $K$ = number of parameters;\n",
                   "\\item AIC = Akaike Information Criterion; BIC = Bayesian Information Criterion.\n",
                   "\\item Pseudo-$R^2$ = McFadden's $R^2 = 1 - LL/LL_0$.\n",
                   "\\end{tablenotes}\n",
                   "\\end{table}\n")

writeLines(tab7_tex, "/Users/hengzhezhao/Desktop/testing/tables/tab_model_comparison.tex")
message("  Saved tab_model_comparison.tex")

message("\n=== All tables generated ===")
