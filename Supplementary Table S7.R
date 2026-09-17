# Supplementary Table S7. Accuracy by correct-answer type.
# Participant-level percent correct; GEE effects are in percentage points.

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

fmt_pct <- function(x) {
  sprintf("%.1f", mean(x, na.rm = TRUE))
}

fmt_pp_ci <- function(est, lo, hi) {
  sprintf("%.2f (%.2f-%.2f)", est, lo, hi)
}

fmt_p <- function(p) {
  if (is.na(p)) return("")
  if (p < 0.001) return("<0.001")
  sprintf("%.3f", p)
}

score_vars <- c(
  "statements1_s",
  "statements2_s",
  "statements3_s",
  "statements4_s",
  "statements5_s",
  "statements6_s",
  "statements7_s",
  "statements8_s",
  "statements9_s",
  "statements10_s"
)

dat <- reviewer_panel(data_path)

missing_vars <- setdiff(
  c(score_vars, "pid", "treatment", "time"),
  names(dat)
)
if (length(missing_vars) > 0) {
  stop("Missing expected columns: ", paste(missing_vars, collapse = ", "))
}

item_values <- unlist(dat[score_vars], use.names = FALSE)
if (!all(item_values %in% c(0, 1))) {
  stop("Item scores are not all 0/1; the correct/incorrect coding changed.")
}

dat <- dat %>%
  mutate(
    treatment = as.integer(treatment),
    time = as.integer(time),
    overall_true = 100 * (statements1_s + statements2_s + statements3_s) / 3,
    overall_false = 100 * (
      statements4_s + statements5_s + statements6_s + statements7_s +
        statements8_s + statements9_s + statements10_s
    ) / 7,
    m1_true = 100 * (statements1_s + statements2_s) / 2,
    m1_false = 100 * (
      statements4_s + statements5_s + statements8_s + statements9_s
    ) / 4,
    m2_true = 100 * statements3_s,
    m2_false = 100 * (
      statements6_s + statements7_s + statements10_s
    ) / 3,
    did = treatment * time,
    id = as.integer(factor(pid))
  ) %>%
  arrange(id, time)

if (!all(dat$treatment %in% 0:1) || !all(dat$time %in% 0:1)) {
  stop("treatment and time must both be coded 0/1")
}

row_specs <- list(
  list(
    module = "Overall",
    type = "TRUE-answer statements",
    var = "overall_true"
  ),
  list(
    module = "",
    type = "FALSE-answer statements",
    var = "overall_false"
  ),
  list(
    module = "Module 1",
    type = "TRUE-answer statements",
    var = "m1_true"
  ),
  list(
    module = "",
    type = "FALSE-answer statements",
    var = "m1_false"
  ),
  list(
    module = "Module 2",
    type = "TRUE-answer statement",
    var = "m2_true"
  ),
  list(
    module = "",
    type = "FALSE-answer statements",
    var = "m2_false"
  )
)

gee_did <- function(df, outcome) {
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
  p <- unname(s["did", "Pr(>|W|)"])
  list(
    estimate = est,
    lo = est - 1.96 * se,
    hi = est + 1.96 * se,
    p = p
  )
}

rows <- lapply(row_specs, function(spec) {
  v <- spec$var
  fc_pre <- dat[[v]][dat$treatment == 1 & dat$time == 0]
  fc_post <- dat[[v]][dat$treatment == 1 & dat$time == 1]
  sv_pre <- dat[[v]][dat$treatment == 0 & dat$time == 0]
  sv_post <- dat[[v]][dat$treatment == 0 & dat$time == 1]
  g <- gee_did(dat, v)

  data.frame(
    Module = spec$module,
    `Statement type` = spec$type,
    `Fact-checking pre, % correct` = fmt_pct(fc_pre),
    `Fact-checking post` = fmt_pct(fc_post),
    `Standard pre` = fmt_pct(sv_pre),
    `Standard post` = fmt_pct(sv_post),
    `Treatment effect, pp (95% CI)` = fmt_pp_ci(g$estimate, g$lo, g$hi),
    P = fmt_p(g$p),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
})

out <- bind_rows(rows)
out_csv <- file.path(out_dir, "Supplementary_Table_S7.csv")
write.csv(out, out_csv, row.names = FALSE)

print(out, row.names = FALSE, right = FALSE)
cat("\nSaved:", out_csv, "\n")
