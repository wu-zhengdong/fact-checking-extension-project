# Supplementary Table S8. Signal-detection sensitivity and response criterion.
# FALSE responses are hits on FALSE-answer items and false alarms on TRUE-answer items.

library(dplyr)
library(tidyr)
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

item_specs <- tribble(
  ~var,           ~module,    ~answer_key,
  "statements2",  "Module 1", "TRUE",
  "statements1",  "Module 1", "TRUE",
  "statements4",  "Module 1", "FALSE",
  "statements8",  "Module 1", "FALSE",
  "statements9",  "Module 1", "FALSE",
  "statements5",  "Module 1", "FALSE",
  "statements6",  "Module 2", "FALSE",
  "statements7",  "Module 2", "FALSE",
  "statements3",  "Module 2", "TRUE",
  "statements10", "Module 2", "FALSE"
)

wide <- reviewer_panel(data_path)

missing_vars <- setdiff(c(item_specs$var, "pid", "treatment", "time"), names(wide))
if (length(missing_vars) > 0) {
  stop("Missing expected columns: ", paste(missing_vars, collapse = ", "))
}

raw_values <- unlist(wide[item_specs$var], use.names = FALSE)
if (!all(raw_values %in% c(1L, 2L, 3L))) {
  stop("Raw statement codes are not all 1/2/3 (True / False / Don't know).")
}

long <- wide %>%
  transmute(
    pid = as.character(pid),
    id = as.integer(factor(pid)),
    treatment = as.integer(treatment),
    time = as.integer(time),
    did = treatment * time,
    across(all_of(item_specs$var), as.integer)
  ) %>%
  pivot_longer(
    cols = all_of(item_specs$var),
    names_to = "item",
    values_to = "raw_response"
  ) %>%
  mutate(
    item = factor(item, levels = item_specs$var),
    answered_false = as.integer(raw_response == 2L)
  ) %>%
  left_join(
    item_specs %>% transmute(item = factor(var, levels = var), module, answer_key),
    by = "item"
  ) %>%
  arrange(id, time, item)

n_participants <- n_distinct(long$id)
stopifnot(nrow(long) == n_participants * 2 * nrow(item_specs))
stopifnot(!anyNA(long$answered_false), !anyNA(long$answer_key))
if (!all(long$treatment %in% 0:1) || !all(long$time %in% 0:1)) {
  stop("treatment and time must both be coded 0/1")
}

# Apply the Hautus correction to every participant-time observation.
# Don't know responses contribute neither a hit nor a false alarm.
hautus_sdt <- function(hits, n_false, false_alarms, n_true) {
  h_rate <- (hits + 0.5) / (n_false + 1)
  fa_rate <- (false_alarms + 0.5) / (n_true + 1)
  z_h <- qnorm(h_rate)
  z_fa <- qnorm(fa_rate)
  tibble(
    n_false = n_false,
    n_true = n_true,
    hits = hits,
    false_alarms = false_alarms,
    hit_rate = h_rate,
    fa_rate = fa_rate,
    dprime = z_h - z_fa,
    criterion = -(z_h + z_fa) / 2
  )
}

sdt_for_items <- function(df, items, scope) {
  sub <- df %>% filter(item %in% items)
  counts <- sub %>%
    group_by(id, treatment, time, did) %>%
    summarise(
      n_false = sum(answer_key == "FALSE"),
      hits = sum(answer_key == "FALSE" & answered_false == 1L),
      n_true = sum(answer_key == "TRUE"),
      false_alarms = sum(answer_key == "TRUE" & answered_false == 1L),
      .groups = "drop"
    )
  expected_n_false <- sum(item_specs$var %in% items & item_specs$answer_key == "FALSE")
  expected_n_true <- sum(item_specs$var %in% items & item_specs$answer_key == "TRUE")
  if (any(counts$n_false != expected_n_false) || any(counts$n_true != expected_n_true)) {
    stop("Item counts for ", scope, " are not constant across participants")
  }
  bind_cols(
    counts %>% dplyr::select(id, treatment, time, did),
    hautus_sdt(counts$hits, counts$n_false, counts$false_alarms, counts$n_true)
  ) %>%
    mutate(scope = scope, .before = 1) %>%
    arrange(id, time)
}

scopes <- list(
  list(name = "Overall", items = item_specs$var),
  list(
    name = "Module 1",
    items = item_specs$var[item_specs$module == "Module 1"]
  ),
  list(
    name = "Module 2",
    items = item_specs$var[item_specs$module == "Module 2"]
  )
)

sdt_list <- lapply(scopes, function(spec) {
  sdt_for_items(long, spec$items, spec$name)
})
sdt_all <- bind_rows(sdt_list)

if (any(!is.finite(sdt_all$dprime)) || any(!is.finite(sdt_all$criterion))) {
  stop("Non-finite d' or c after Hautus correction")
}
if (nrow(sdt_all) != n_participants * 2 * length(scopes)) {
  stop("Unexpected number of participant-time-module rows")
}

fmt_rate <- function(x) sprintf("%.3f", mean(x, na.rm = TRUE))

fmt_mean <- function(x) sprintf("%.2f", mean(x, na.rm = TRUE))

fmt_pp_ci <- function(est, lo, hi) {
  sprintf("%.2f (%.2f-%.2f)", est, lo, hi)
}

fmt_p <- function(p) {
  if (is.na(p)) return("")
  if (p < 0.001) return("<0.001")
  sprintf("%.3f", p)
}

fmt <- function(x, digits = 4) formatC(x, format = "f", digits = digits)

cell <- function(df, outcome, treatment_value, time_value, formatter) {
  x <- df[[outcome]][df$treatment == treatment_value & df$time == time_value]
  formatter(x)
}

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
    se = se,
    lo = est - 1.96 * se,
    hi = est + 1.96 * se,
    p = p
  )
}

make_scope_rows <- function(df, scope) {
  df <- df %>% filter(scope == !!scope) %>% arrange(id, time)
  g_d <- gee_did(df, "dprime")
  g_c <- gee_did(df, "criterion")

  measures <- list(
    list(
      measure = "Hit rate",
      outcome = "hit_rate",
      formatter = fmt_rate,
      effect = "—",
      p = "—"
    ),
    list(
      measure = "False-alarm rate",
      outcome = "fa_rate",
      formatter = fmt_rate,
      effect = "—",
      p = "—"
    ),
    list(
      measure = "d'",
      outcome = "dprime",
      formatter = fmt_mean,
      effect = fmt_pp_ci(g_d$estimate, g_d$lo, g_d$hi),
      p = fmt_p(g_d$p)
    ),
    list(
      measure = "Criterion c",
      outcome = "criterion",
      formatter = fmt_mean,
      effect = fmt_pp_ci(g_c$estimate, g_c$lo, g_c$hi),
      p = fmt_p(g_c$p)
    )
  )

  rows <- bind_rows(lapply(seq_along(measures), function(i) {
    spec <- measures[[i]]
    tibble(
      Module = if (i == 1L) scope else "",
      Measure = spec$measure,
      `Fact-checking pre` = cell(df, spec$outcome, 1L, 0L, spec$formatter),
      `Fact-checking post` = cell(df, spec$outcome, 1L, 1L, spec$formatter),
      `Standard-view pre` = cell(df, spec$outcome, 0L, 0L, spec$formatter),
      `Standard-view post` = cell(df, spec$outcome, 0L, 1L, spec$formatter),
      `Treatment effect (95% CI)` = spec$effect,
      P = spec$p
    )
  }))

  details <- tibble(
    Module = scope,
    Section = c(
      rep("Items", 2),
      rep("Hautus correction", 2),
      rep("DiD d'", 5),
      rep("DiD criterion c", 5),
      "Sample"
    ),
    Quantity = c(
      "FALSE-answer items (hits)",
      "TRUE-answer items (false alarms)",
      "H",
      "FA",
      "Estimate", "Standard error", "95% CI lower", "95% CI upper", "P value",
      "Estimate", "Standard error", "95% CI lower", "95% CI upper", "P value",
      "Participant-time rows"
    ),
    Value = c(
      as.character(df$n_false[1]),
      as.character(df$n_true[1]),
      "(hits + 0.5) / (n_false + 1)",
      "(false_alarms + 0.5) / (n_true + 1)",
      fmt(g_d$estimate), fmt(g_d$se), fmt(g_d$lo), fmt(g_d$hi),
      if (g_d$p < 0.0001) "<0.0001" else fmt(g_d$p),
      fmt(g_c$estimate), fmt(g_c$se), fmt(g_c$lo), fmt(g_c$hi),
      if (g_c$p < 0.0001) "<0.0001" else fmt(g_c$p),
      format(nrow(df), big.mark = ",")
    )
  )

  list(table = rows, details = details)
}

blocks <- lapply(c("Overall", "Module 1", "Module 2"), function(scope) {
  make_scope_rows(sdt_all, scope)
})

out <- bind_rows(lapply(blocks, function(b) b$table))
details <- bind_rows(lapply(blocks, function(b) b$details))

out_csv <- file.path(out_dir, "Supplementary_Table_S8.csv")
out_details <- file.path(out_dir, "details/Supplementary_Table_S8_details.csv")
write.csv(out, file.path(out_dir, "details", "Supplementary_Table_S8_full.csv"), row.names = FALSE)
for (i in seq_len(nrow(out))) if (out$Module[i] == "") out$Module[i] <- out$Module[i - 1L]
out <- out[out$Measure %in% c("d'", "Criterion c"), ]
out <- out[order(match(out$Module, c("Module 1", "Module 2", "Overall"))), ]
write.csv(out, out_csv, row.names = FALSE)
write.csv(details, out_details, row.names = FALSE)

print(as.data.frame(out), row.names = FALSE, right = FALSE)
cat("\nSaved:", out_csv, "\n")
cat("Saved:", out_details, "\n")
