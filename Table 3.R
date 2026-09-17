# Table 3. Vaccine confidence and vaccination willingness.
# Positive response: 4 or 5 on the five-point scale.

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
dir.create(file.path(script_dir, "analysis_tables", "details/Table_3"), recursive = TRUE, showWarnings = FALSE)

fmt_n_pct <- function(n, denom) {
  sprintf("%d (%.1f)", n, if (denom > 0) 100 * n / denom else 0)
}

fmt_rr_ci <- function(est, lo, hi) {
  sprintf("%.2f (%.2f-%.2f)", est, lo, hi)
}

fmt_p <- function(p) {
  if (is.na(p)) return("")
  if (p < 0.001) return("<0.001")
  sprintf("%.3f", p)
}

chi_p_time <- function(df, var, arm) {
  sub <- df %>% filter(treatment == arm)
  ct <- table(sub$time, sub[[var]])
  if (nrow(ct) < 2 || ncol(ct) < 2) return(NA_real_)
  suppressWarnings(chisq.test(ct, correct = FALSE)$p.value)
}

gee_rr_did <- function(df, outcome) {
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
  p <- unname(s["did", "Pr(>|W|)"])
  list(
    rr = exp(b),
    lo = exp(b - 1.96 * se),
    hi = exp(b + 1.96 * se),
    p = p
  )
}

dat <- read_xlsx(data_path, sheet = "panel") %>%
  mutate(
    conf_1_1_new = as.integer(conf_1_new >= 4),
    conf_2_1_new = as.integer(conf_2_new >= 4),
    willingness_1_new = as.integer(willingness_new >= 4),
    did = treatment * time,
    id = as.integer(factor(pid))
  ) %>%
  arrange(id, time)

stopifnot(!anyNA(dat$conf_1_1_new), !anyNA(dat$conf_2_1_new), !anyNA(dat$willingness_1_new))

row_specs <- list(
  list(
    label = "Perceived vaccine efficacy (positive)",
    var = "conf_1_1_new",
    note = "conf_1_new = 4 or 5"
  ),
  list(
    label = "Perceived vaccine safety (positive)",
    var = "conf_2_1_new",
    note = "conf_2_new = 4 or 5"
  ),
  list(
    label = "Willingness to vaccinate daughter (positive)",
    var = "willingness_1_new",
    note = "willingness_new = 4 or 5"
  )
)

rows <- lapply(row_specs, function(spec) {
  v <- spec$var
  n_fc0 <- sum(dat$treatment == 1 & dat$time == 0)
  n_fc1 <- sum(dat$treatment == 1 & dat$time == 1)
  n_sv0 <- sum(dat$treatment == 0 & dat$time == 0)
  n_sv1 <- sum(dat$treatment == 0 & dat$time == 1)

  k_fc0 <- sum(dat[[v]] == 1 & dat$treatment == 1 & dat$time == 0)
  k_fc1 <- sum(dat[[v]] == 1 & dat$treatment == 1 & dat$time == 1)
  k_sv0 <- sum(dat[[v]] == 1 & dat$treatment == 0 & dat$time == 0)
  k_sv1 <- sum(dat[[v]] == 1 & dat$treatment == 0 & dat$time == 1)

  p_fc <- chi_p_time(dat, v, arm = 1)
  p_sv <- chi_p_time(dat, v, arm = 0)
  g <- gee_rr_did(dat, v)

  data.frame(
    Outcome = spec$label,
    `Fact-checking Pre, n (%)` = fmt_n_pct(k_fc0, n_fc0),
    `Fact-checking Post, n (%)` = fmt_n_pct(k_fc1, n_fc1),
    `Fact-checking P (chi2)` = fmt_p(p_fc),
    `Standard view Pre, n (%)` = fmt_n_pct(k_sv0, n_sv0),
    `Standard view Post, n (%)` = fmt_n_pct(k_sv1, n_sv1),
    `Standard view P (chi2)` = fmt_p(p_sv),
    `RR (95% CI)` = fmt_rr_ci(g$rr, g$lo, g$hi),
    `P value` = fmt_p(g$p),
    Coding = spec$note,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
})

table3 <- bind_rows(rows)

out_csv <- file.path(out_dir, "Table3.csv")
write.csv(table3, out_csv, row.names = FALSE)

cat("Saved:", out_csv, "\n")
