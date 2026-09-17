# Supplementary Table S4. Item-level logistic mixed-effects models.
# Random intercepts for participants and items; treatment-by-time effect reported as an OR.

library(dplyr)
library(tidyr)
library(lme4)

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

item_specs <- tribble(
  ~var,             ~module,    ~section,
  "statements2_s",  "Module 1", "Correct information exposure",
  "statements1_s",  "Module 1", "Correct information exposure",
  "statements4_s",  "Module 1", "Correct information exposure",
  "statements8_s",  "Module 1", "Misinformation exposure",
  "statements9_s",  "Module 1", "Misinformation exposure",
  "statements5_s",  "Module 1", "Misinformation exposure",
  "statements6_s",  "Module 2", "Correct information exposure",
  "statements7_s",  "Module 2", "Correct information exposure",
  "statements3_s",  "Module 2", "Misinformation exposure",
  "statements10_s", "Module 2", "Misinformation exposure"
)

wide <- reviewer_panel(data_path)

missing_vars <- setdiff(c(item_specs$var, "pid", "treatment", "time"), names(wide))
if (length(missing_vars) > 0) {
  stop("Missing expected columns: ", paste(missing_vars, collapse = ", "))
}

item_values <- unlist(wide[item_specs$var], use.names = FALSE)
if (!all(item_values %in% c(0, 1))) {
  stop("Item scores are not all 0/1; the correct/incorrect coding changed.")
}

long <- wide %>%
  transmute(
    pid = factor(pid),
    treatment = as.integer(treatment),
    time = as.integer(time),
    across(all_of(item_specs$var), as.integer)
  ) %>%
  pivot_longer(
    cols = all_of(item_specs$var),
    names_to = "item",
    values_to = "correct"
  ) %>%
  mutate(
    item = factor(item, levels = item_specs$var),
    correct = as.integer(correct)
  ) %>%
  arrange(pid, time, item)

n_participants <- nlevels(long$pid)
n_items <- nlevels(long$item)
n_obs <- nrow(long)

rows_per_id <- table(long$pid)
if (any(rows_per_id != n_items * 2)) {
  stop(
    "Expected ", n_items * 2, " item responses per participant; ",
    sum(rows_per_id != n_items * 2), " participants differ."
  )
}
stopifnot(n_obs == n_participants * n_items * 2, !anyNA(long$correct))
if (!all(long$treatment %in% 0:1) || !all(long$time %in% 0:1)) {
  stop("treatment and time must both be coded 0/1")
}

score_check <- long %>%
  group_by(pid, time) %>%
  summarise(score_from_long = sum(correct), .groups = "drop") %>%
  inner_join(
    wide %>%
      transmute(
        pid = factor(pid),
        time = as.integer(time),
        score_wide = rowSums(across(all_of(item_specs$var)))
      ),
    by = c("pid", "time")
  )
stopifnot(
  nrow(score_check) == n_participants * 2,
  all(score_check$score_from_long == score_check$score_wide)
)

fmt_pct <- function(data, treatment_value, time_value) {
  x <- data$correct[data$treatment == treatment_value & data$time == time_value]
  sprintf("%.1f", 100 * mean(x))
}

fmt_or_ci <- function(est, lo, hi) {
  sprintf("%.2f (%.2f-%.2f)", est, lo, hi)
}

INTERACTION_TERM <- "treatment:time"

fmt_p <- function(p) {
  if (is.na(p)) return("")
  if (p < 0.001) return("<0.001")
  sprintf("%.3f", p)
}

fmt <- function(x, digits = 4) formatC(x, format = "f", digits = digits)

fit_item_model <- function(data, label) {
  data <- droplevels(data)

  fit <- glmer(
    correct ~ treatment * time + (1 | pid) + (1 | item),
    data = data,
    family = binomial(link = "logit")
  )

  coefs <- summary(fit)$coefficients
  if (!INTERACTION_TERM %in% rownames(coefs)) {
    stop("Interaction term ", INTERACTION_TERM, " not found for ", label)
  }

  est <- unname(coefs[INTERACTION_TERM, "Estimate"])
  se <- unname(coefs[INTERACTION_TERM, "Std. Error"])
  p_value <- unname(coefs[INTERACTION_TERM, "Pr(>|z|)"])
  ci <- suppressMessages(
    confint(fit, parm = INTERACTION_TERM, method = "Wald")
  )
  ci_low <- unname(ci[1, 1])
  ci_high <- unname(ci[1, 2])

  var_components <- as.data.frame(VarCorr(fit))
  var_pid <- var_components$vcov[var_components$grp == "pid"]
  var_item <- var_components$vcov[var_components$grp == "item"]

  opt_info <- fit@optinfo
  conv_code <- opt_info$conv$opt
  conv_messages <- opt_info$conv$lme4$messages
  n_item <- nlevels(data$item)

  list(
    fit = fit,
    label = label,
    n_items = n_item,
    fc_pre = fmt_pct(data, 1L, 0L),
    fc_post = fmt_pct(data, 1L, 1L),
    sv_pre = fmt_pct(data, 0L, 0L),
    sv_post = fmt_pct(data, 0L, 1L),
    or = exp(est),
    or_lo = exp(ci_low),
    or_hi = exp(ci_high),
    p_value = p_value,
    details = tibble(
      Outcome = label,
      Section = c(
        rep("Model", 3),
        rep("Fixed effect (treatment x time)", 8),
        rep("Random effects", 4),
        rep("Convergence", 5),
        rep("Sample", 3)
      ),
      Quantity = c(
        "Formula", "Family", "Estimation",
        "Term", "Coefficient (log-odds)", "Standard error", "Odds ratio",
        "OR 95% CI lower", "OR 95% CI upper", "95% CI method", "P value",
        "Participant (pid) intercept variance", "Participant (pid) intercept SD",
        "Item (item) intercept variance", "Item (item) intercept SD",
        "Optimizer", "Convergence code (0 = converged)", "Converged",
        "Convergence warnings", "Singular fit",
        "Participants", "Items", "Item-response observations"
      ),
      Value = c(
        "correct ~ treatment * time + (1 | pid) + (1 | item)",
        "binomial(link = \"logit\")",
        "lme4::glmer, Laplace approximation (nAGQ = 1)",
        INTERACTION_TERM, fmt(est), fmt(se), fmt(exp(est)),
        fmt(exp(ci_low)), fmt(exp(ci_high)), "Wald",
        if (p_value < 0.0001) "<0.0001" else fmt(p_value),
        fmt(var_pid), fmt(sqrt(var_pid)), fmt(var_item), fmt(sqrt(var_item)),
        paste(opt_info$optimizer, collapse = ", "),
        as.character(conv_code),
        if (isTRUE(conv_code == 0) && length(conv_messages) == 0) "Yes" else "No",
        if (length(conv_messages) == 0) "None" else paste(conv_messages, collapse = "; "),
        if (isSingular(fit)) "Yes" else "No",
        format(nlevels(data$pid), big.mark = ","),
        as.character(n_item),
        format(nrow(data), big.mark = ",")
      )
    )
  )
}

model_specs <- list(
  list(label = "Overall", items = item_specs$var),
  list(label = "Module 1", items = item_specs$var[item_specs$module == "Module 1"]),
  list(label = "Module 2", items = item_specs$var[item_specs$module == "Module 2"])
)

models <- lapply(model_specs, function(spec) {
  fit_item_model(long %>% filter(item %in% spec$items), spec$label)
})

summary_table <- bind_rows(lapply(models, function(m) {
  tibble(
    Analysis = sprintf("%s, %d items", m$label, m$n_items),
    `Fact-checking pre, % correct` = m$fc_pre,
    `Fact-checking post, % correct` = m$fc_post,
    `Standard-view pre, % correct` = m$sv_pre,
    `Standard-view post, % correct` = m$sv_post,
    `OR (95% CI)` = fmt_or_ci(m$or, m$or_lo, m$or_hi),
    `P value` = fmt_p(m$p_value)
  )
}))

details_table <- bind_rows(lapply(models, function(m) m$details))

out_summary <- file.path(out_dir, "Supplementary_Table_S4.csv")
out_details <- file.path(out_dir, "details/Supplementary_Table_S4_details.csv")
write.csv(summary_table, out_summary, row.names = FALSE)
write.csv(details_table, out_details, row.names = FALSE)

cat("\nTable 2 item-response-level logistic mixed-effects sensitivity analysis\n")
cat(strrep("=", 70), "\n")
for (m in models) {
  cat("\n### ", m$label, " (", m$n_items, " items)\n\n", sep = "")
  print(summary(m$fit))
}
cat("\n", strrep("=", 70), "\n\nSummary\n\n", sep = "")
print(as.data.frame(summary_table), row.names = FALSE, right = FALSE)
cat("\nSaved:", out_summary, "\n")
cat("Saved:", out_details, "\n")
