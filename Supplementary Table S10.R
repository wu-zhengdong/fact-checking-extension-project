# Supplementary Table S10. Class-clustered sensitivity analysis of attitude outcomes.

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
dir.create(file.path(script_dir, "analysis_tables", "details/Supplementary_Table_S10"), recursive = TRUE, showWarnings = FALSE)

fmt_rr <- function(x) sprintf("%.2f", x)

fmt_ci <- function(lo, hi) sprintf("%.2f-%.2f", lo, hi)

fmt_n_pct <- function(n, denom) {
  sprintf("%d (%.1f)", n, if (denom > 0) 100 * n / denom else 0)
}

fmt_p <- function(p) {
  if (is.na(p)) return("")
  if (p < 0.001) return("<0.001")
  sprintf("%.3f", p)
}

dat <- read_xlsx(data_path, sheet = "panel") %>%
  mutate(
    conf_1_1_new = as.integer(conf_1_new >= 4),
    conf_2_1_new = as.integer(conf_2_new >= 4),
    willingness_1_new = as.integer(willingness_new >= 4),
    did = treatment * time,
    pid_int = as.integer(factor(pid)),
    class_int = as.integer(factor(class_id))
  )

stopifnot(
  !anyNA(dat$conf_1_1_new),
  !anyNA(dat$conf_2_1_new),
  !anyNA(dat$willingness_1_new)
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
cat("Table 3 sensitivity analysis: GEE clustered on class instead of pid\n")
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
  list(
    label = "Perceived vaccine efficacy (positive)",
    paper_label = "HPV vaccines are effective, n (%)",
    var = "conf_1_1_new",
    note = "conf_1_new = 4 or 5"
  ),
  list(
    label = "Perceived vaccine safety (positive)",
    paper_label = "HPV vaccines are safe, n (%)",
    var = "conf_2_1_new",
    note = "conf_2_new = 4 or 5"
  ),
  list(
    label = "Willingness to vaccinate daughter (positive)",
    paper_label = "Parental willingness to vaccinate their daughter, n (%)",
    var = "willingness_1_new",
    note = "willingness_new = 4 or 5"
  )
)

gee_rr_did <- function(df, outcome, corstr) {
  fit <- tryCatch(
    geeglm(
      as.formula(paste(outcome, "~ treatment + time + did")),
      id = cluster_id,
      data = df,
      family = binomial(link = "log"),
      corstr = corstr
    ),
    error = function(e) {
      warning(
        "GEE failed for ", outcome, " (", corstr, "): ",
        conditionMessage(e),
        call. = FALSE
      )
      NULL
    }
  )
  if (is.null(fit)) {
    return(list(
      rr = NA_real_, lo = NA_real_, hi = NA_real_,
      se = NA_real_, p = NA_real_, log_rr = NA_real_
    ))
  }
  s <- summary(fit)$coefficients
  b <- unname(s["did", "Estimate"])
  se <- unname(s["did", "Std.err"])
  list(
    log_rr = b,
    rr = exp(b),
    lo = exp(b - 1.96 * se),
    hi = exp(b + 1.96 * se),
    se = se,
    p = unname(s["did", "Pr(>|W|)"])
  )
}

rows <- lapply(row_specs, function(spec) {
  v <- spec$var
  g_pid <- gee_rr_did(dat_pid, v, "exchangeable")
  g_class <- gee_rr_did(dat_class, v, "independence")
  g_class_exch <- gee_rr_did(dat_class, v, "exchangeable")

  corstr_gap <- if (is.na(g_class$log_rr) || is.na(g_class_exch$log_rr)) {
    NA_real_
  } else {
    max(
      abs(g_class$log_rr - g_class_exch$log_rr),
      abs(g_class$se - g_class_exch$se)
    )
  }

  se_ratio <- if (is.na(g_pid$se) || is.na(g_class$se) || g_pid$se == 0) {
    NA_real_
  } else {
    g_class$se / g_pid$se
  }

  cat(sprintf(
    "%-22s RR %.4f | SE(log) pid %.4f -> class %.4f (x%.2f) | exch gap %.1e\n",
    v, g_class$rr, g_pid$se, g_class$se, se_ratio, corstr_gap
  ))

  n_fc0 <- sum(dat$treatment == 1 & dat$time == 0)
  n_fc1 <- sum(dat$treatment == 1 & dat$time == 1)
  n_sv0 <- sum(dat$treatment == 0 & dat$time == 0)
  n_sv1 <- sum(dat$treatment == 0 & dat$time == 1)

  data.frame(
    Outcome = spec$label,
    `Paper label` = spec$paper_label,
    `Fact-checking Pre, n (%)` = fmt_n_pct(
      sum(dat[[v]] == 1 & dat$treatment == 1 & dat$time == 0), n_fc0
    ),
    `Fact-checking Post, n (%)` = fmt_n_pct(
      sum(dat[[v]] == 1 & dat$treatment == 1 & dat$time == 1), n_fc1
    ),
    `Standard view Pre, n (%)` = fmt_n_pct(
      sum(dat[[v]] == 1 & dat$treatment == 0 & dat$time == 0), n_sv0
    ),
    `Standard view Post, n (%)` = fmt_n_pct(
      sum(dat[[v]] == 1 & dat$treatment == 0 & dat$time == 1), n_sv1
    ),
    `Original GEE RR` = fmt_rr(g_pid$rr),
    `Original GEE 95% CI` = fmt_ci(g_pid$lo, g_pid$hi),
    `Original GEE P value` = fmt_p(g_pid$p),
    `Class-clustered GEE RR` = fmt_rr(g_class$rr),
    `Class-clustered GEE 95% CI` = fmt_ci(g_class$lo, g_class$hi),
    `Class-clustered GEE P value` = fmt_p(g_class$p),
    `Number of classes` = n_classes,
    `Number of participants` = n_participants,
    `Fact-checking group n` = n_fc0,
    `Standard view group n` = n_sv0,
    Coding = spec$note,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
})

supp_table3 <- bind_rows(rows)
supp_csv <- file.path(out_dir, "Supplementary_Table_S10.csv")
write.csv(supp_table3, supp_csv, row.names = FALSE)

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
icc_csv <- file.path(out_dir, "details/Supplementary_Table_S10/class_level_icc_table3.csv")
write.csv(icc, icc_csv, row.names = FALSE)

cat("\n")
print(icc[, c("Outcome", "ICC")], row.names = FALSE)
cat("\nSaved:", supp_csv, "\n")
cat("Saved:", icc_csv, "\n")
cat("Main analysis results/table3_attitudes.csv left untouched.\n")
