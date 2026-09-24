# Table 2. Information-assessment scores.
# Participant-clustered linear GEE with exchangeable working correlation.
# Standardized effects use the pooled baseline SD, as in the manuscript.

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

fmt_mean_sd <- function(x) {
  x <- x[!is.na(x)]
  sprintf("%.2f (%.1f)", mean(x), stats::sd(x))
}

fmt_coef_ci <- function(est, lo, hi) {
  sprintf("%.2f (%.2f-%.2f)", est, lo, hi)
}

fmt_p <- function(p) {
  if (is.na(p)) return("")
  if (p < 0.001) return("<0.001")
  sprintf("%.3f", p)
}

items <- read_reviewer_items(data_path)

stopifnot(nrow(items) == 1548L * 2L * 10L)
stopifnot(n_distinct(items$participant_id) == 1548L)
if (!all(items$response %in% c("TRUE", "FALSE", "Don't know"))) {
  stop("response is not the original three-category coding")
}
derived_accuracy <- as.integer(items$response == items$correct_answer)
if (!isTRUE(all(derived_accuracy == items$accuracy))) {
  stop("accuracy is not response == correct_answer (Don't know = 0)")
}

dat <- items %>%
  mutate(
    treatment = as.integer(group == "fact-checking"),
    time_num = as.integer(time == "post")
  ) %>%
  group_by(participant_id, group, treatment, time, time_num) %>%
  summarise(
    score1 = sum(accuracy[analysis_module == "Module 1" & post_veracity == "accurate"]),
    score2 = sum(accuracy[analysis_module == "Module 1" & post_veracity == "misinformation"]),
    score4 = sum(accuracy[analysis_module == "Module 2" & post_veracity == "accurate"]),
    score5 = sum(accuracy[analysis_module == "Module 2" & post_veracity == "misinformation"]),
    n_items = n(),
    .groups = "drop"
  ) %>%
  mutate(
    score3 = score1 + score2,
    score6 = score4 + score5,
    score7 = score3 + score6,
    time = time_num,
    did = treatment * time,
    id = as.integer(factor(participant_id))
  ) %>%
  arrange(id, time)

stopifnot(nrow(dat) == 1548L * 2L, all(dat$n_items == 10L))
stopifnot(all(dat$score1 %in% 0:3), all(dat$score2 %in% 0:3))
stopifnot(all(dat$score4 %in% 0:2), all(dat$score5 %in% 0:2))

row_specs <- list(
  list(module = "Module 1", exposure = "Correct information", var = "score1"),
  list(module = "", exposure = "Misinformation", var = "score2"),
  list(module = "", exposure = "Score (0-6)", var = "score3"),
  list(module = "Module 2", exposure = "Correct information", var = "score4"),
  list(module = "", exposure = "Misinformation", var = "score5"),
  list(module = "", exposure = "Score (0-4)", var = "score6"),
  list(module = "Overall score (0-10)", exposure = "", var = "score7")
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
  g <- gee_did(dat, v)
  baseline_1 <- dat[[v]][dat$treatment == 1 & dat$time == 0]
  baseline_0 <- dat[[v]][dat$treatment == 0 & dat$time == 0]
  n1 <- length(baseline_1)
  n0 <- length(baseline_0)
  pooled_sd <- sqrt(((n1 - 1) * var(baseline_1) +
                     (n0 - 1) * var(baseline_0)) / (n1 + n0 - 2))
  if (!is.finite(pooled_sd) || pooled_sd <= 0) {
    stop("Cannot standardize ", v, ": pooled baseline SD must be positive.")
  }
  # Treat the baseline pooled SD as a fixed standardiser: rescale the
  # GEE robust confidence limits, matching table2_effect_sizes.R.
  cohens_d <- fmt_coef_ci(g$estimate / pooled_sd,
                         g$lo / pooled_sd, g$hi / pooled_sd)
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
    `Coefficient (95% CI)` = fmt_coef_ci(g$estimate, g$lo, g$hi),
    `P value` = fmt_p(g$p),
    `Cohen's d (95% CI)` = cohens_d,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
})

table2 <- bind_rows(rows)
out_csv <- file.path(out_dir, "Table2.csv")
write.csv(table2, out_csv, row.names = FALSE)

print(table2, row.names = FALSE, right = FALSE)
cat("Saved:", out_csv, "\n")
