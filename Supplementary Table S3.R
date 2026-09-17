# Supplementary Table S3. Item-specific accuracy in Module 2.

library(dplyr)
library(geepack)

args_cmd <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args_cmd, value = TRUE)
script_dir <- if (length(file_arg)) dirname(normalizePath(gsub("~+~", " ", sub("^--file=", "", file_arg[1]), fixed = TRUE))) else {
  getwd()
}
source(file.path(script_dir, "lib", "reviewer_data.R"))
root <- script_dir
out_dir <- file.path(root, "analysis_tables")
dir.create(file.path(out_dir, "details"), recursive = TRUE, showWarnings = FALSE)
data_path <- reviewer_data_path(root)

dat <- reviewer_panel(data_path) %>%
  mutate(
    treatment = as.integer(treatment),
    time = as.integer(time),
    did = treatment * time,
    id = as.integer(factor(pid))
  ) %>%
  arrange(id, time)

fmt_n_pct <- function(n, denom) {
  sprintf("%d (%.1f)", n, 100 * n / denom)
}

fmt_p <- function(p) {
  if (is.na(p)) return("")
  if (p < 0.001) return("<0.001")
  sprintf("%.3f", p)
}

fmt_rr_ci <- function(rr, lo, hi) {
  sprintf("%.2f (%.2f-%.2f)", rr, lo, hi)
}

chi_p_time <- function(df, outcome, group_var, group_value) {
  sub <- df[df[[group_var]] == group_value, , drop = FALSE]
  ct <- table(sub$time, sub[[outcome]])
  if (nrow(ct) < 2 || ncol(ct) < 2) return(NA_real_)
  suppressWarnings(chisq.test(ct, correct = FALSE)$p.value)
}

gee_rr <- function(df, outcome, group_var, did_var) {
  model_formula <- as.formula(
    paste(outcome, "~", group_var, "+ time +", did_var)
  )
  fit <- tryCatch(
    geeglm(
      model_formula,
      id = id,
      data = df,
      family = binomial(link = "log"),
      corstr = "exchangeable"
    ),
    error = function(e) NULL
  )
  # Use modified Poisson if the log-binomial fit fails.
  if (is.null(fit)) {
    message("Log-binomial model failed for ", outcome, "; using modified Poisson")
    fit <- geeglm(
      model_formula,
      id = id,
      data = df,
      family = poisson(link = "log"),
      corstr = "exchangeable"
    )
  }
  s <- summary(fit)$coefficients
  b <- unname(s[did_var, "Estimate"])
  se <- unname(s[did_var, "Std.err"])
  p <- unname(s[did_var, "Pr(>|W|)"])
  list(
    rr = exp(b),
    lo = exp(b - 1.96 * se),
    hi = exp(b + 1.96 * se),
    p = p
  )
}

make_item_table <- function(df, specs, group_var, did_var, group_labels) {
  g1 <- 1L
  g0 <- 0L
  n1 <- sum(df[[group_var]] == g1 & df$time == 0)
  n0 <- sum(df[[group_var]] == g0 & df$time == 0)

  rows <- lapply(specs, function(spec) {
    outcome <- spec$var
    k_1_pre <- sum(df[[outcome]] == 1 & df[[group_var]] == g1 & df$time == 0)
    k_1_post <- sum(df[[outcome]] == 1 & df[[group_var]] == g1 & df$time == 1)
    k_0_pre <- sum(df[[outcome]] == 1 & df[[group_var]] == g0 & df$time == 0)
    k_0_post <- sum(df[[outcome]] == 1 & df[[group_var]] == g0 & df$time == 1)

    p1 <- chi_p_time(df, outcome, group_var, g1)
    p0 <- chi_p_time(df, outcome, group_var, g0)
    g <- gee_rr(df, outcome, group_var, did_var)

    row <- data.frame(
      `Exposure information` = spec$section,
      Statement = spec$label,
      group1_pre = fmt_n_pct(k_1_pre, n1),
      group1_post = fmt_n_pct(k_1_post, n1),
      group1_p = fmt_p(p1),
      group0_pre = fmt_n_pct(k_0_pre, n0),
      group0_post = fmt_n_pct(k_0_post, n0),
      group0_p = fmt_p(p0),
      `RR (95% CI)` = fmt_rr_ci(g$rr, g$lo, g$hi),
      `P value` = fmt_p(g$p),
      `Group 1 N` = n1,
      `Group 0 N` = n0,
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
    names(row)[3:8] <- c(
      paste0(group_labels[["one"]], " Pre, n (%)"),
      paste0(group_labels[["one"]], " Post, n (%)"),
      paste0(group_labels[["one"]], " P value"),
      paste0(group_labels[["zero"]], " Pre, n (%)"),
      paste0(group_labels[["zero"]], " Post, n (%)"),
      paste0(group_labels[["zero"]], " P value")
    )
    row
  })

  bind_rows(rows)
}

module2_specs <- list(
  list(
    section = "Correct information exposure",
    var = "statements6_s",
    label = paste(
      "Women who are breastfeeding should not get the HPV vaccine because",
      "it can make the milk have vaccine components, which is not good",
      "for the baby."
    )
  ),
  list(
    section = "Correct information exposure",
    var = "statements7_s",
    label = "Cervical cancer screening no longer needed after HPV vaccination."
  ),
  list(
    section = "Misinformation exposure",
    var = "statements3_s",
    label = "Bivalent HPV vaccine is enough, early vaccination is more important."
  ),
  list(
    section = "Misinformation exposure",
    var = "statements10_s",
    label = "HPV vaccination often leads to severe allergic reactions."
  )
)

module2 <- make_item_table(
  dat,
  module2_specs,
  group_var = "treatment",
  did_var = "did",
  group_labels = c(
    one = "Fact-checking extension",
    zero = "Standard view"
  )
)

write.csv(module2, file.path(out_dir, "Supplementary_Table_S3.csv"), row.names = FALSE)
