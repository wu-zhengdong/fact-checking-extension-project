# Table 1. Baseline sample characteristics.

library(readxl)
library(dplyr)

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
dir.create(file.path(script_dir, "analysis_tables", "details/Table_1"), recursive = TRUE, showWarnings = FALSE)

region_levels <- c(3, 1, 2) # Megacity, Urban, Rural (paper order)
region_labels <- c(
  "3" = "Megacity",
  "1" = "Urban",
  "2" = "Rural"
)

grade_labels <- c(
  "1" = "Grade 6",
  "2" = "Grade 7",
  "3" = "Grade 8",
  "4" = "Grade 9"
)

fmt_n_pct <- function(n, denom) {
  sprintf("%d (%.1f)", n, if (denom > 0) 100 * n / denom else 0)
}

fmt_mean_sd <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) return("")
  sprintf("%.1f (%.1f)", mean(x), sd(x))
}

fmt_p <- function(p) {
  if (is.na(p)) return("")
  if (p < 0.001) return("<0.001")
  sprintf("%.3f", p)
}

empty_row <- function(characteristic = "", total = "", fc = "", sv = "", p = "") {
  data.frame(
    Characteristics = characteristic,
    Total = total,
    `Fact-checking` = fc,
    `Standard view` = sv,
    `P value` = p,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

dat <- read_xlsx(data_path, sheet = "panel") %>%
  filter(time == 0) %>%
  mutate(
    treatment = as.integer(treatment),
    district_1 = as.integer(district_1),
    grade = as.integer(grade)
  )

n_total <- nrow(dat)
n_fc <- sum(dat$treatment == 1)
n_sv <- sum(dat$treatment == 0)

stopifnot(n_total == 1548L, n_fc == 771L, n_sv == 777L)

chi_p <- function(df, var) {
  ct <- table(df[[var]], df$treatment)
  if (nrow(ct) < 2 || ncol(ct) < 2) return(NA_real_)
  suppressWarnings(chisq.test(ct, correct = FALSE)$p.value)
}

ttest_p <- function(df, var) {
  x1 <- df[[var]][df$treatment == 1]
  x0 <- df[[var]][df$treatment == 0]
  x1 <- x1[!is.na(x1)]
  x0 <- x0[!is.na(x0)]
  if (length(x1) < 2 || length(x0) < 2) return(NA_real_)
  t.test(x1, x0, var.equal = TRUE)$p.value
}

add_categorical <- function(df, var, label, levels, level_labels) {
  p <- chi_p(df, var)
  rows <- list(empty_row(sprintf("%s, n (%%)", label), p = fmt_p(p)))

  for (lv in levels) {
    lab <- unname(level_labels[as.character(lv)])
    if (is.na(lab) || is.null(lab)) lab <- as.character(lv)

    n_t <- sum(df[[var]] == lv, na.rm = TRUE)
    n_fc_i <- sum(df[[var]] == lv & df$treatment == 1, na.rm = TRUE)
    n_sv_i <- sum(df[[var]] == lv & df$treatment == 0, na.rm = TRUE)

    rows[[length(rows) + 1]] <- empty_row(
      characteristic = paste0("  ", lab),
      total = fmt_n_pct(n_t, n_total),
      fc = fmt_n_pct(n_fc_i, n_fc),
      sv = fmt_n_pct(n_sv_i, n_sv)
    )
  }
  bind_rows(rows)
}

add_continuous <- function(df, var, label) {
  p <- ttest_p(df, var)
  empty_row(
    characteristic = sprintf("%s, mean (SD)", label),
    total = fmt_mean_sd(df[[var]]),
    fc = fmt_mean_sd(df[[var]][df$treatment == 1]),
    sv = fmt_mean_sd(df[[var]][df$treatment == 0]),
    p = fmt_p(p)
  )
}

parts <- list(
  empty_row("Participants' daughter"),

  add_categorical(
    dat, "district_1", "Region",
    levels = region_levels,
    level_labels = region_labels
  ),

  add_categorical(
    dat, "grade", "Grade",
    levels = c(1, 2, 3, 4),
    level_labels = grade_labels
  ),

  add_continuous(dat, "age", "Age (years)"),

  empty_row("Participants"),

  add_categorical(
    dat, "gender_parents", "Parenthood",
    levels = c(0, 1),
    level_labels = c("0" = "Mother", "1" = "Father")
  ),

  add_continuous(dat, "age_parents", "Parent age (years)"),

  add_categorical(
    dat, "education_2", "Education",
    levels = c(1, 2, 3),
    level_labels = c(
      "1" = "Junior middle school and below",
      "2" = "Senior high school",
      "3" = "College and above"
    )
  ),

  add_categorical(
    dat, "employment_1", "Formal employment",
    levels = c(1, 0),
    level_labels = c("1" = "Yes", "0" = "No")
  ),

  add_categorical(
    dat, "income_1", "Annual household income",
    levels = c(1, 2, 3, 4),
    level_labels = c(
      "1" = "<100k CNY",
      "2" = "100-200k CNY",
      "3" = "200-300k CNY",
      "4" = ">300k CNY"
    )
  ),

  add_categorical(
    dat, "maternal_vaccination", "Self/spouse HPV vaccination",
    levels = c(1, 0),
    level_labels = c("1" = "Yes", "0" = "No")
  ),

  add_categorical(
    dat, "fata_exp_1", "Past-year HPV misinformation exposure",
    levels = c(0, 1),
    level_labels = c("0" = "No", "1" = "Yes")
  )
)

table1 <- bind_rows(parts)
colnames(table1) <- c(
  "Characteristics",
  sprintf("Total (N = %s)", format(n_total, big.mark = ",")),
  sprintf("Fact-checking (N = %s)", format(n_fc, big.mark = ",")),
  sprintf("Standard view (N = %s)", format(n_sv, big.mark = ",")),
  "P value"
)

out_csv <- file.path(out_dir, "Table1.csv")

write.csv(table1, out_csv, row.names = FALSE)

cat("Saved:", out_csv, "\n")
cat(sprintf("N: Total=%d, Fact-checking=%d, Standard view=%d\n", n_total, n_fc, n_sv))
