# Table 5. User experience among participants exposed to the extension.
# Includes all 771 fact-checking participants and 537 standard-view participants who opted in.

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
dir.create(file.path(script_dir, "analysis_tables", "details/Table_5"), recursive = TRUE, showWarnings = FALSE)

col_outcome <- "User-experience outcome"

col_header <- function(label, n) {
  sprintf("%s (N = %s), n (%%)", label, format(n, big.mark = ","))
}

fmt_n_pct <- function(n, denom) {
  sprintf("%d (%.1f)", n, 100 * n / denom)
}

fmt_mean_sd <- function(x) {
  sprintf("%.2f (%.2f)", mean(x), sd(x))
}

fmt_median_iqr <- function(x) {
  q <- quantile(x, c(0.25, 0.5, 0.75), type = 7)
  sprintf("%.0f (%.0f-%.0f)", q[2], q[1], q[3])
}

items <- tibble::tribble(
  ~label,                                                         ~var,
  "Satisfied with the extension",                                 "exten_satisfication",
  "Perceived that the extension enhanced HPV vaccine confidence",  "exten_confidence",
  "Willing to use the extension",                                  "exten_willingness",
  "Willing to recommend the extension",                            "exten_recommen"
)

dat <- read_xlsx(data_path, sheet = "panel") %>%
  filter(time == 1) %>%
  mutate(
    exposure = case_when(
      treatment == 1 ~ "fact_checking",
      treatment == 0 & switch_3_1 == 1 ~ "standard_view_opt_in",
      TRUE ~ "standard_view_declined"
    )
  )

stopifnot(
  nrow(dat) == 1548L,
  sum(dat$treatment == 1) == 771L,
  sum(dat$treatment == 0) == 777L,
  !any(is.na(dat$switch_3_1)),
  !any(is.na(dat[items$var]))
)

exposed <- dat %>% filter(exposure != "standard_view_declined")

n_all <- nrow(exposed)
n_fc <- sum(exposed$exposure == "fact_checking")
n_sv <- sum(exposed$exposure == "standard_view_opt_in")
n_excluded <- sum(dat$exposure == "standard_view_declined")

stopifnot(n_fc == 771L, n_sv == 537L, n_all == 1308L, n_excluded == 240L)

col_all <- col_header("All exposed participants", n_all)
col_fc <- col_header("Fact-checking group", n_fc)
col_sv <- col_header("Standard-view opt-in participants", n_sv)

is_positive <- function(x) x %in% c(4, 5)

subsets <- list(
  all = exposed,
  fc = exposed %>% filter(exposure == "fact_checking"),
  sv = exposed %>% filter(exposure == "standard_view_opt_in")
)

table4_desc <- purrr::pmap_dfr(items, function(label, var) {
  cell <- function(d) fmt_n_pct(sum(is_positive(d[[var]])), nrow(d))
  out <- data.frame(
    label,
    cell(subsets$all),
    cell(subsets$fc),
    cell(subsets$sv),
    stringsAsFactors = FALSE
  )
  names(out) <- c(col_outcome, col_all, col_fc, col_sv)
  out
})

out_csv <- file.path(out_dir, "Table5.csv")
write.csv(table4_desc, out_csv, row.names = FALSE)

group_labels <- c(
  all = "All exposed participants",
  fc = "Fact-checking group",
  sv = "Standard-view opt-in participants"
)

details <- purrr::pmap_dfr(items, function(label, var) {
  purrr::imap_dfr(subsets, function(d, key) {
    x <- d[[var]]
    data.frame(
      Outcome = label,
      Variable = var,
      Group = group_labels[[key]],
      Definition = c(
        all = "treatment == 1 | (treatment == 0 & switch_3_1 == 1)",
        fc = "treatment == 1",
        sv = "treatment == 0 & switch_3_1 == 1"
      )[[key]],
      N = nrow(d),
      `Positive (score 4 or 5), n` = sum(is_positive(x)),
      `Positive (score 4 or 5), %` = sprintf("%.1f", 100 * mean(is_positive(x))),
      `Score 1, n` = sum(x == 1),
      `Score 2, n` = sum(x == 2),
      `Score 3, n` = sum(x == 3),
      `Score 4, n` = sum(x == 4),
      `Score 5, n` = sum(x == 5),
      `Mean (SD)` = fmt_mean_sd(x),
      `Median (IQR)` = fmt_median_iqr(x),
      Missing = sum(is.na(x)),
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  })
})

out_details <- file.path(out_dir, "details/Table_5/table4_user_experience_descriptive_details.csv")
write.csv(details, out_details, row.names = FALSE)

enabled <- dat %>% filter(switch_3_1 == 1)
enabled_subsets <- list(
  all = enabled,
  fc = enabled %>% filter(treatment == 1),
  sv = enabled %>% filter(treatment == 0)
)
stopifnot(
  nrow(enabled_subsets$all) == 1100L,
  nrow(enabled_subsets$fc) == 563L,
  nrow(enabled_subsets$sv) == 537L
)

sens <- purrr::pmap_dfr(items, function(label, var) {
  cell <- function(d) fmt_n_pct(sum(is_positive(d[[var]])), nrow(d))
  out <- data.frame(
    label,
    cell(enabled_subsets$all),
    cell(enabled_subsets$fc),
    cell(enabled_subsets$sv),
    stringsAsFactors = FALSE
  )
  names(out) <- c(
    col_outcome,
    col_header("Enabled the extension", nrow(enabled_subsets$all)),
    col_header("Fact-checking group who enabled", nrow(enabled_subsets$fc)),
    col_header("Standard-view opt-in participants", nrow(enabled_subsets$sv))
  )
  out
})

out_sens <- file.path(
  out_dir, "details/Table_5/table4_user_experience_descriptive_sensitivity.csv"
)
write.csv(sens, out_sens, row.names = FALSE)

print_md <- function(df) {
  cells <- lapply(df, as.character)
  widths <- mapply(
    function(x, nm) max(nchar(c(nm, x))),
    cells, names(df)
  )
  pad <- function(x, w, left) formatC(x, width = if (left) -w else w)
  cat("|", paste(
    mapply(pad, names(df), widths, c(TRUE, rep(FALSE, ncol(df) - 1))),
    collapse = " | "
  ), "|\n", sep = " ")
  cat("|", paste(mapply(function(w, left) {
    if (left) paste0(strrep("-", w + 1)) else paste0(strrep("-", w), ":")
  }, widths, c(TRUE, rep(FALSE, ncol(df) - 1))), collapse = "|"), "|\n", sep = "")
  for (i in seq_len(nrow(df))) {
    row <- mapply(
      function(x, w, left) pad(x[i], w, left),
      cells, widths, c(TRUE, rep(FALSE, ncol(df) - 1))
    )
    cat("|", paste(row, collapse = " | "), "|\n", sep = " ")
  }
}

cat("\nTable 4 (descriptive). Exposure strata, post-exposure survey:\n")
cat(sprintf("  All exposed participants          N = %s\n", format(n_all, big.mark = ",")))
cat(sprintf("  Fact-checking group               N = %s\n", format(n_fc, big.mark = ",")))
cat(sprintf("  Standard-view opt-in participants N = %s\n", format(n_sv, big.mark = ",")))
cat(sprintf(
  "  Excluded (standard view, extension never enabled) N = %s\n",
  format(n_excluded, big.mark = ",")
))
cat("\nPositive response = score 4 or 5 on the 1-5 agreement scale.\n\n")
print_md(table4_desc)

cat("\nSensitivity: columns restricted to participants who enabled the",
    "extension at Q30/Q43\n")
cat(sprintf(
  "  N = %s enabled (%s fact-checking, %s standard view)\n\n",
  format(nrow(enabled_subsets$all), big.mark = ","),
  format(nrow(enabled_subsets$fc), big.mark = ","),
  format(nrow(enabled_subsets$sv), big.mark = ",")
))
print_md(sens)

cat("\nSaved:", out_csv, "\n")
cat("Saved:", out_details, "\n")
cat("Saved:", out_sens, "\n")
