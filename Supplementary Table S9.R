# Supplementary Table S9. Sensitivity analysis using the five-point attitude scales.
# Linear GEE and proportional-odds models; higher scores indicate more positive responses.

library(readxl)
library(dplyr)
library(tidyr)
library(geepack)
if (!requireNamespace("MASS", quietly = TRUE)) stop("Package 'MASS' is required.")

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
dir.create(file.path(script_dir, "analysis_tables", "details/Supplementary_Table_S9"), recursive = TRUE, showWarnings = FALSE)

LIKERT_LEVELS <- 1:5
# Resample participants, retaining both observations in each bootstrap sample.
N_BOOT <- 2000
BOOT_SEED <- 20260827

fmt <- function(x, digits = 4) formatC(x, format = "f", digits = digits)

fmt_mean_sd <- function(x) sprintf("%.2f (%.2f)", mean(x), stats::sd(x))

fmt_median_iqr <- function(x) {
  q <- stats::quantile(x, c(0.25, 0.5, 0.75), type = 7)
  sprintf("%.0f [%.0f-%.0f]", q[2], q[1], q[3])
}

fmt_est_ci <- function(est, lo, hi) sprintf("%.2f (%.2f to %.2f)", est, lo, hi)

fmt_ratio_ci <- function(est, lo, hi) sprintf("%.2f (%.2f-%.2f)", est, lo, hi)

fmt_p <- function(p) {
  if (is.na(p)) return("")
  if (p < 0.001) return("<0.001")
  sprintf("%.3f", p)
}

outcome_specs <- list(
  list(
    label = "Perceived vaccine efficacy",
    var = "conf_1_new",
    binary_var = "conf_1_1_new",
    item = "Q13 (baseline) / Q50 (follow-up), stored as 6 - raw",
    baseline_src = "conf_1",
    followup_src = "conf_1_follow",
    low = "1 = not effective at all",
    high = "5 = very effective"
  ),
  list(
    label = "Perceived vaccine safety",
    var = "conf_2_new",
    binary_var = "conf_2_1_new",
    item = "Q14 (baseline) / Q51 (follow-up), stored as 6 - raw",
    baseline_src = "conf_2",
    followup_src = "conf_2_follow",
    low = "1 = not safe at all",
    high = "5 = very safe"
  ),
  list(
    label = "Willingness to vaccinate daughter",
    var = "willingness_new",
    binary_var = "willingness_1_new",
    item = "Q16 (baseline) / Q52 (follow-up), stored as 6 - raw",
    baseline_src = "willingness",
    followup_src = "willingness_follow",
    low = "1 = definitely unwilling",
    high = "5 = definitely willing"
  )
)

likert_vars <- vapply(outcome_specs, function(s) s$var, character(1))

raw <- read_xlsx(data_path, sheet = "panel")

missing_vars <- setdiff(c(likert_vars, "pid", "treatment", "time"), names(raw))
if (length(missing_vars) > 0) {
  stop("Missing expected columns: ", paste(missing_vars, collapse = ", "))
}

dat <- raw %>%
  mutate(
    conf_1_1_new = as.integer(conf_1_new >= 4),
    conf_2_1_new = as.integer(conf_2_new >= 4),
    willingness_1_new = as.integer(willingness_new >= 4),
    did = treatment * time,
    id = as.integer(factor(pid))
  ) %>%
  arrange(id, time)

for (v in likert_vars) {
  values <- dat[[v]]
  if (anyNA(values)) stop("Missing values in ", v)
  if (!all(values %in% LIKERT_LEVELS)) {
    stop(v, " has values outside 1-5; the response coding changed.")
  }
}

stopifnot(
  all(dat$conf_1_1_new == as.integer(dat$conf_1_new >= 4)),
  all(dat$conf_2_1_new == as.integer(dat$conf_2_new >= 4)),
  all(dat$willingness_1_new == as.integer(dat$willingness_new >= 4))
)

for (spec in outcome_specs) {
  src_cols <- c(spec$baseline_src, spec$followup_src)
  if (!all(src_cols %in% names(dat))) {
    stop("Cannot verify scale direction for ", spec$var,
         ": missing ", paste(setdiff(src_cols, names(dat)), collapse = ", "))
  }
  at_baseline <- dat$time == 0
  if (!isTRUE(all.equal(dat[[spec$var]][at_baseline], dat[[spec$baseline_src]][at_baseline])) ||
      !isTRUE(all.equal(dat[[spec$var]][!at_baseline], dat[[spec$followup_src]][!at_baseline]))) {
    stop(
      spec$var, " does not match ", spec$baseline_src, "/", spec$followup_src,
      "; the scale direction assumed by this script can no longer be verified."
    )
  }
}

n_participants <- length(unique(dat$id))
rows_per_id <- table(dat$id)
if (any(rows_per_id != 2)) {
  stop(
    "Expected 2 observations per participant; ",
    sum(rows_per_id != 2), " participants differ."
  )
}
stopifnot(nrow(dat) == n_participants * 2, all(dat$did == dat$treatment * dat$time))

gee_linear_did <- function(df, outcome) {
  fit <- geeglm(
    as.formula(paste(outcome, "~ treatment + time + did")),
    id = id,
    data = df,
    family = gaussian(link = "identity"),
    corstr = "exchangeable"
  )
  s <- summary(fit)$coefficients
  est <- unname(s["did", "Estimate"])
  se <- unname(s["did", "Std.err"])
  list(
    fit = fit,
    estimate = est,
    se = se,
    lo = est - 1.96 * se,
    hi = est + 1.96 * se,
    p = unname(s["did", "Pr(>|W|)"]),
    alpha = unname(summary(fit)$corr["alpha", "Estimate"])
  )
}

gee_log_binomial_did <- function(df, outcome) {
  fit <- tryCatch(
    geeglm(
      as.formula(paste(outcome, "~ treatment + time + did")),
      id = id,
      data = df,
      family = binomial(link = "log"),
      corstr = "exchangeable"
    ),
    error = function(e) {
      warning("GEE failed for ", outcome, ": ", conditionMessage(e), call. = FALSE)
      NULL
    }
  )
  if (is.null(fit)) {
    return(list(rr = NA_real_, lo = NA_real_, hi = NA_real_, p = NA_real_))
  }
  s <- summary(fit)$coefficients
  b <- unname(s["did", "Estimate"])
  se <- unname(s["did", "Std.err"])
  list(
    rr = exp(b),
    lo = exp(b - 1.96 * se),
    hi = exp(b + 1.96 * se),
    p = unname(s["did", "Pr(>|W|)"])
  )
}

polr_did_coef <- function(df, outcome) {
  model_df <- data.frame(
    y = factor(df[[outcome]], levels = LIKERT_LEVELS, ordered = TRUE),
    treatment = df$treatment,
    time = df$time,
    did = df$did
  )
  model_df$y <- droplevels(model_df$y)
  if (nlevels(model_df$y) < 2) return(NA_real_)
  fit <- tryCatch(
    MASS::polr(y ~ treatment + time + did, data = model_df, Hess = FALSE, method = "logistic"),
    error = function(e) NULL
  )
  if (is.null(fit) || !("did" %in% names(coef(fit)))) return(NA_real_)
  unname(coef(fit)["did"])
}

polr_cluster_bootstrap <- function(df, outcome, n_boot, seed) {
  observed <- polr_did_coef(df, outcome)
  if (is.na(observed)) {
    stop("Proportional-odds model did not fit for ", outcome)
  }

  ids <- unique(df$id)
  rows_by_id <- split(seq_len(nrow(df)), df$id)

  set.seed(seed)
  draws <- vapply(seq_len(n_boot), function(i) {
    sampled <- sample(ids, size = length(ids), replace = TRUE)
    boot_df <- df[unlist(rows_by_id[as.character(sampled)], use.names = FALSE), , drop = FALSE]
    polr_did_coef(boot_df, outcome)
  }, numeric(1))

  n_failed <- sum(is.na(draws))
  draws <- draws[!is.na(draws)]
  if (length(draws) < 0.9 * n_boot) {
    stop(
      "Only ", length(draws), " of ", n_boot,
      " bootstrap resamples converged for ", outcome
    )
  }

  se <- stats::sd(draws)
  ci <- unname(stats::quantile(draws, c(0.025, 0.975), type = 7))

  list(
    estimate = observed,
    se = se,
    or = exp(observed),
    or_lo = exp(ci[1]),
    or_hi = exp(ci[2]),
    p = 2 * stats::pnorm(-abs(observed / se)),
    n_boot = n_boot,
    n_used = length(draws),
    n_failed = n_failed
  )
}

arm_slice <- function(df, outcome, arm, tp) {
  df[[outcome]][df$treatment == arm & df$time == tp]
}

within_arm_tests <- function(df, outcome, arm) {
  pre <- arm_slice(df, outcome, arm, 0)
  post <- arm_slice(df, outcome, arm, 1)
  list(
    t = stats::t.test(post, pre, var.equal = TRUE)$p.value,
    wilcox = suppressWarnings(stats::wilcox.test(post, pre)$p.value)
  )
}

distribution_table <- bind_rows(lapply(outcome_specs, function(spec) {
  dat %>%
    mutate(
      Arm = factor(
        ifelse(treatment == 1, "Fact-checking", "Standard view"),
        levels = c("Fact-checking", "Standard view")
      ),
      Period = factor(ifelse(time == 1, "Post", "Pre"), levels = c("Pre", "Post")),
      response = factor(.data[[spec$var]], levels = LIKERT_LEVELS)
    ) %>%
    count(Arm, Period, response, .drop = FALSE) %>%
    group_by(Arm, Period) %>%
    mutate(pct = 100 * n / sum(n)) %>%
    ungroup() %>%
    arrange(Arm, Period) %>%
    transmute(
      Outcome = spec$label,
      Variable = spec$var,
      Arm,
      Period,
      Response = as.integer(as.character(response)),
      `n (%)` = sprintf("%d (%.1f)", n, pct)
    )
})) %>%
  pivot_wider(
    names_from = Response,
    values_from = `n (%)`,
    names_prefix = "Response "
  ) %>%
  rename(
    `Response 1 (least positive), n (%)` = `Response 1`,
    `Response 2, n (%)` = `Response 2`,
    `Response 3, n (%)` = `Response 3`,
    `Response 4, n (%)` = `Response 4`,
    `Response 5 (most positive), n (%)` = `Response 5`
  )

cat("Fitting Likert-scale sensitivity models (", N_BOOT,
    " bootstrap resamples per outcome)\n", sep = "")

results <- lapply(outcome_specs, function(spec) {
  cat("  ", spec$label, "...\n", sep = "")
  v <- spec$var

  linear <- gee_linear_did(dat, v)
  ordinal <- polr_cluster_bootstrap(dat, v, n_boot = N_BOOT, seed = BOOT_SEED)
  binary <- gee_log_binomial_did(dat, spec$binary_var)

  fc <- within_arm_tests(dat, v, arm = 1)
  sv <- within_arm_tests(dat, v, arm = 0)

  list(
    spec = spec,
    linear = linear,
    ordinal = ordinal,
    binary = binary,
    fc = fc,
    sv = sv,
    fc_pre = arm_slice(dat, v, 1, 0),
    fc_post = arm_slice(dat, v, 1, 1),
    sv_pre = arm_slice(dat, v, 0, 0),
    sv_post = arm_slice(dat, v, 0, 1)
  )
})

summary_table <- bind_rows(lapply(results, function(r) {
  data.frame(
    Outcome = r$spec$label,
    `Fact-checking Pre, mean (s.d.)` = fmt_mean_sd(r$fc_pre),
    `Fact-checking Post, mean (s.d.)` = fmt_mean_sd(r$fc_post),
    `Fact-checking P (t-test)` = fmt_p(r$fc$t),
    `Standard view Pre, mean (s.d.)` = fmt_mean_sd(r$sv_pre),
    `Standard view Post, mean (s.d.)` = fmt_mean_sd(r$sv_post),
    `Standard view P (t-test)` = fmt_p(r$sv$t),
    `Mean difference (95% CI)` = fmt_est_ci(r$linear$estimate, r$linear$lo, r$linear$hi),
    `P value (linear GEE)` = fmt_p(r$linear$p),
    `OR (95% CI)` = fmt_ratio_ci(r$ordinal$or, r$ordinal$or_lo, r$ordinal$or_hi),
    `P value (ordinal)` = fmt_p(r$ordinal$p),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}))

comparison_table <- bind_rows(lapply(results, function(r) {
  data.frame(
    Outcome = r$spec$label,
    `Main analysis: RR (95% CI), positive response 4-5` =
      fmt_ratio_ci(r$binary$rr, r$binary$lo, r$binary$hi),
    `Main analysis: P value` = fmt_p(r$binary$p),
    `Sensitivity: mean difference (95% CI), 1-5 scale` =
      fmt_est_ci(r$linear$estimate, r$linear$lo, r$linear$hi),
    `Sensitivity: P value` = fmt_p(r$linear$p),
    `Sensitivity: OR (95% CI), proportional odds` =
      fmt_ratio_ci(r$ordinal$or, r$ordinal$or_lo, r$ordinal$or_hi),
    `Sensitivity: P value (ordinal)` = fmt_p(r$ordinal$p),
    `Same direction` = ifelse(
      is.na(r$binary$rr), "",
      ifelse(sign(r$binary$rr - 1) == sign(r$linear$estimate) &&
               sign(r$linear$estimate) == sign(r$ordinal$estimate), "Yes", "No")
    ),
    `Same conclusion at 0.05` = ifelse(
      is.na(r$binary$p), "",
      ifelse((r$binary$p < 0.05) == (r$linear$p < 0.05) &&
               (r$linear$p < 0.05) == (r$ordinal$p < 0.05), "Yes", "No")
    ),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}))

details_table <- bind_rows(lapply(results, function(r) {
  lin <- r$linear
  ord <- r$ordinal
  data.frame(
    Outcome = r$spec$label,
    Section = c(
      rep("Outcome", 6),
      rep("Linear GEE (1-5 scale)", 9),
      rep("Proportional odds (1-5 scale)", 10),
      rep("Descriptives", 6),
      rep("Within-arm tests", 4),
      rep("Sample", 2)
    ),
    Quantity = c(
      "Variable", "Source item", "Response scale",
      "Lowest category", "Highest category", "Dichotomized counterpart in table3.R",
      "Formula", "Family", "Working correlation", "Estimated alpha",
      "Term", "Mean difference (scale points; + = more positive)",
      "Standard error", "95% CI", "P value",
      "Formula", "Method", "Inference",
      "Term", "Coefficient (log-odds)", "Bootstrap standard error",
      "Odds ratio (OR > 1 = shift to more positive)",
      "OR 95% CI (percentile)", "P value",
      "Bootstrap resamples converged",
      "Fact-checking pre, mean (s.d.)", "Fact-checking post, mean (s.d.)",
      "Fact-checking post, median [IQR]",
      "Standard view pre, mean (s.d.)", "Standard view post, mean (s.d.)",
      "Standard view post, median [IQR]",
      "Fact-checking P (two-sample t-test)", "Fact-checking P (Wilcoxon rank-sum)",
      "Standard view P (two-sample t-test)", "Standard view P (Wilcoxon rank-sum)",
      "Participants", "Observations"
    ),
    Value = c(
      r$spec$var, r$spec$item,
      "1-5 Likert (raw, not dichotomized); higher = more positive",
      r$spec$low, r$spec$high, r$spec$binary_var,
      paste(r$spec$var, "~ treatment + time + did"),
      "gaussian(link = \"identity\")", "exchangeable", fmt(lin$alpha),
      "did", fmt(lin$estimate), fmt(lin$se),
      sprintf("%s to %s", fmt(lin$lo), fmt(lin$hi)),
      if (lin$p < 0.0001) "<0.0001" else fmt(lin$p),
      sprintf("ordered(%s) ~ treatment + time + did", r$spec$var),
      "MASS::polr, logistic (proportional odds)",
      sprintf("Nonparametric bootstrap, %d resamples of participants, seed %d",
              ord$n_boot, BOOT_SEED),
      "did", fmt(ord$estimate), fmt(ord$se),
      fmt(ord$or), sprintf("%s to %s", fmt(ord$or_lo), fmt(ord$or_hi)),
      if (ord$p < 0.0001) "<0.0001" else fmt(ord$p),
      sprintf("%d of %d", ord$n_used, ord$n_boot),
      fmt_mean_sd(r$fc_pre), fmt_mean_sd(r$fc_post), fmt_median_iqr(r$fc_post),
      fmt_mean_sd(r$sv_pre), fmt_mean_sd(r$sv_post), fmt_median_iqr(r$sv_post),
      fmt(r$fc$t), fmt(r$fc$wilcox), fmt(r$sv$t), fmt(r$sv$wilcox),
      format(n_participants, big.mark = ","), format(nrow(dat), big.mark = ",")
    ),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}))

out_summary <- file.path(out_dir, "Supplementary_Table_S9.csv")
out_comparison <- file.path(out_dir, "details/Supplementary_Table_S9/table3_likert_sensitivity_comparison.csv")
out_distribution <- file.path(out_dir, "details/Supplementary_Table_S9/table3_likert_sensitivity_distribution.csv")
out_details <- file.path(out_dir, "details/Supplementary_Table_S9/table3_likert_sensitivity_details.csv")

write.csv(summary_table, out_summary, row.names = FALSE)
write.csv(comparison_table, out_comparison, row.names = FALSE)
write.csv(distribution_table, out_distribution, row.names = FALSE)
write.csv(details_table, out_details, row.names = FALSE)

cat("\nLikert-scale sensitivity analysis for Table 3\n")
cat(strrep("=", 70), "\n")
for (r in results) {
  cat("\n### ", r$spec$label, " (", r$spec$var, ")\n\n", sep = "")
  print(summary(r$linear$fit))
}
cat("\n", strrep("=", 70), "\n\n", sep = "")
cat("Scale direction: 1 = least positive, 5 = most positive ",
    "(the analysis file stores 6 - raw item).\n", sep = "")
cat("A positive mean difference, or an OR above 1, means the fact-checking arm ",
    "moved toward the positive end.\n\n", sep = "")
cat("Summary (raw 1-5 scale)\n\n")
print(as.data.frame(summary_table), row.names = FALSE, right = FALSE)
cat("\nDichotomized main analysis versus Likert sensitivity analysis\n\n")
print(as.data.frame(comparison_table), row.names = FALSE, right = FALSE)
cat("\nSaved:", out_summary, "\n")
cat("Saved:", out_comparison, "\n")
cat("Saved:", out_distribution, "\n")
cat("Saved:", out_details, "\n")
