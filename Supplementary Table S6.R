# Supplementary Table S6. Class-clustered sensitivity analysis of information-assessment scores.

library(readxl)
library(dplyr)
library(geepack)

args_cmd <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args_cmd, value = TRUE)
script_dir <- if (length(file_arg)) dirname(normalizePath(gsub("~+~", " ", sub("^--file=", "", file_arg[1]), fixed = TRUE))) else {
  getwd()
}
root <- file.path(script_dir, "data-code-table")
project_root <- root

data_path <- file.path(root, "analysis_data", "processed", "data_full.xlsx")
out_dir <- file.path(script_dir, "analysis_tables")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(script_dir, "analysis_tables", "details/Supplementary_Table_S6"), recursive = TRUE, showWarnings = FALSE)

fmt_coef <- function(x) sprintf("%.2f", x)

fmt_ci <- function(lo, hi) sprintf("%.2f-%.2f", lo, hi)

fmt_mean_sd <- function(x) {
  x <- x[!is.na(x)]
  sprintf("%.2f (%.1f)", mean(x), stats::sd(x))
}

fmt_p <- function(p) {
  if (is.na(p)) return("")
  if (p < 0.001) return("<0.001")
  sprintf("%.3f", p)
}

dat <- read_xlsx(data_path, sheet = "panel") %>%
  mutate(
    score1 = statements2_s + statements1_s + statements4_s,
    score2 = statements8_s + statements9_s + statements5_s,
    score3 = score1 + score2,
    score4 = statements6_s + statements7_s,
    score5 = statements3_s + statements10_s,
    score6 = score4 + score5,
    score7 = score3 + score6,
    did = treatment * time,
    pid_int = as.integer(factor(pid)),
    class_int = as.integer(factor(class_id))
  )

if (!"class_id" %in% names(dat)) {
  stop("class_id not found; rebuild data_full.xlsx with export_data_full.py")
}
if (any(is.na(dat$class_id))) stop("Missing class_id values")

n_classes <- n_distinct(dat$class_id)
n_participants <- n_distinct(dat$pid)

pid_per_class <- dat %>%
  distinct(pid, class_id) %>%
  count(pid, name = "n_classes")
if (any(pid_per_class$n_classes != 1)) {
  stop("A participant maps to more than one class_id")
}

class_sizes <- dat %>%
  filter(time == 0) %>%
  count(class_id, name = "n")

cat(strrep("=", 70), "\n", sep = "")
cat("Table 2 sensitivity analysis: GEE clustered on class instead of pid\n")
cat(strrep("=", 70), "\n", sep = "")
cat("Participants:", n_participants, "| classes:", n_classes, "\n")
cat(sprintf(
  "Baseline participants per class - min %d, median %.1f, mean %.2f, max %d\n",
  min(class_sizes$n), stats::median(class_sizes$n),
  mean(class_sizes$n), max(class_sizes$n)
))

# GEE requires contiguous rows within each cluster.
dat_pid <- dat %>%
  arrange(pid_int, time) %>%
  mutate(cluster_id = pid_int)
dat_class <- dat %>%
  arrange(class_int, pid_int, time) %>%
  mutate(cluster_id = class_int)

row_specs <- list(
  list(module = "Module 1", exposure = "Correct information", var = "score1"),
  list(module = "", exposure = "Misinformation", var = "score2"),
  list(module = "", exposure = "Score (1-6)", var = "score3"),
  list(module = "Module 2", exposure = "Correct information", var = "score4"),
  list(module = "", exposure = "Misinformation", var = "score5"),
  list(module = "", exposure = "Score (1-4)", var = "score6"),
  list(module = "Overall score (1-10)", exposure = "", var = "score7")
)

gee_did <- function(df, outcome, corstr) {
  fit <- geeglm(
    as.formula(paste(outcome, "~ treatment + time + did")),
    id = cluster_id,
    data = df,
    family = gaussian(link = "identity"),
    corstr = corstr
  )
  s <- summary(fit)$coefficients
  est <- unname(s["did", "Estimate"])
  se <- unname(s["did", "Std.err"])
  list(
    estimate = est,
    lo = est - 1.96 * se,
    hi = est + 1.96 * se,
    se = se,
    p = unname(s["did", "Pr(>|W|)"])
  )
}

CORSTR_TOL <- 1e-8

rows <- lapply(row_specs, function(spec) {
  v <- spec$var
  g_pid <- gee_did(dat_pid, v, "exchangeable")
  g_class <- gee_did(dat_class, v, "independence")

  g_class_exch <- gee_did(dat_class, v, "exchangeable")
  corstr_gap <- max(
    abs(g_class$estimate - g_class_exch$estimate),
    abs(g_class$se - g_class_exch$se)
  )
  if (corstr_gap > CORSTR_TOL) {
    stop(sprintf(
      "%s: independence and exchangeable disagree by %.2e at the class level",
      v, corstr_gap
    ))
  }

  cat(sprintf(
    "%-8s DiD %.4f | SE pid %.4f -> class %.4f (x%.2f) | exch gap %.1e\n",
    v, g_class$estimate, g_pid$se, g_class$se, g_class$se / g_pid$se,
    corstr_gap
  ))

  data.frame(
    Module = spec$module,
    `Exposure information` = spec$exposure,
    `Fact-checking Pre-exposure score, mean (s.d.)` =
      fmt_mean_sd(dat[[v]][dat$treatment == 1 & dat$time == 0]),
    `Fact-checking Post-exposure score, mean (s.d.)` =
      fmt_mean_sd(dat[[v]][dat$treatment == 1 & dat$time == 1]),
    `Standard view Pre-exposure score, mean (s.d.)` =
      fmt_mean_sd(dat[[v]][dat$treatment == 0 & dat$time == 0]),
    `Standard view Post-exposure score, mean (s.d.)` =
      fmt_mean_sd(dat[[v]][dat$treatment == 0 & dat$time == 1]),
    `Original GEE estimate` = fmt_coef(g_pid$estimate),
    `Original GEE 95% CI` = fmt_ci(g_pid$lo, g_pid$hi),
    `Original GEE P value` = fmt_p(g_pid$p),
    `Class-clustered GEE estimate` = fmt_coef(g_class$estimate),
    `Class-clustered GEE 95% CI` = fmt_ci(g_class$lo, g_class$hi),
    `Class-clustered GEE P value` = fmt_p(g_class$p),
    `Number of classes` = n_classes,
    `Number of participants` = n_participants,
    `Fact-checking group n` = sum(dat$treatment == 1 & dat$time == 0),
    `Standard view group n` = sum(dat$treatment == 0 & dat$time == 0),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
})

supp_table2 <- bind_rows(rows)
supp_csv <- file.path(out_dir, "Supplementary_Table_S6.csv")
write.csv(supp_table2, supp_csv, row.names = FALSE)

# Baseline ICC from a random-intercept model.
baseline <- dat %>% filter(time == 0)
use_lme4 <- requireNamespace("lme4", quietly = TRUE)
cat("\nICC fitted with:", if (use_lme4) "lme4::lmer" else "nlme::lme", "\n")

icc_one <- function(outcome) {
  form <- as.formula(paste(outcome, "~ 1"))
  if (use_lme4) {
    fit <- lme4::lmer(
      as.formula(paste(outcome, "~ 1 + (1 | class_id)")),
      data = baseline, REML = TRUE
    )
    vc <- as.data.frame(lme4::VarCorr(fit))
    var_class <- vc$vcov[vc$grp == "class_id"]
    var_resid <- vc$vcov[vc$grp == "Residual"]
  } else {
    fit <- nlme::lme(
      form, random = ~ 1 | class_id, data = baseline, method = "REML"
    )
    sds <- as.numeric(nlme::VarCorr(fit)[, "StdDev"])
    var_class <- sds[1]^2
    var_resid <- sds[2]^2
  }
  data.frame(
    Outcome = outcome,
    `Class variance` = sprintf("%.4f", var_class),
    `Residual variance` = sprintf("%.4f", var_resid),
    ICC = sprintf("%.4f", var_class / (var_class + var_resid)),
    `Number of classes` = n_classes,
    `Number of participants` = n_participants,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

icc <- bind_rows(lapply(vapply(row_specs, function(s) s$var, ""), icc_one))
icc_csv <- file.path(out_dir, "details/Supplementary_Table_S6/class_level_icc.csv")
write.csv(icc, icc_csv, row.names = FALSE)

cat("\n")
print(icc[, c("Outcome", "ICC")], row.names = FALSE)
cat("\nSaved:", supp_csv, "\n")
cat("Saved:", icc_csv, "\n")
cat("Main analysis results/table2_module_scores.csv left untouched.\n")
