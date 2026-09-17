# Table 6. Subgroup analyses with multiplicity-adjusted interaction P values.

library(geepack)
library(dplyr)
library(purrr)
library(readxl)
library(stringr)
library(tibble)

# Adjust eight subgroup tests within each outcome ("outcome"), or all 24 ("all").
adjust_method <- "BH"
adjust_family <- "outcome"

format_est_ci <- function(est, low, high, digits = 2) {
  sprintf(
    paste0("%.", digits, "f (%.", digits, "f, %.", digits, "f)"),
    est, low, high
  )
}

format_p <- function(p) {
  if (is.na(p)) return("")
  format.pval(p, digits = 3, eps = 0.001)
}

safe_solve <- function(vcov_matrix) {
  tryCatch(
    solve(vcov_matrix),
    error = function(e) MASS::ginv(vcov_matrix)
  )
}

get_arm_time_effect <- function(data, outcome) {
  empty <- list(estimate = "", p_value = NA_real_)

  data_model <- data %>%
    filter(
      !is.na(.data[[outcome]]),
      !is.na(arm),
      !is.na(time)
    ) %>%
    droplevels()

  if (n_distinct(data_model$arm) < 2 || n_distinct(data_model$time) < 2) {
    return(empty)
  }

  fit <- tryCatch(
    geeglm(
      as.formula(paste0(outcome, " ~ arm * time")),
      id = id,
      data = data_model,
      family = gaussian(link = "identity"),
      corstr = "exchangeable"
    ),
    error = function(e) {
      warning("GEE failed for ", outcome, ": ", conditionMessage(e), call. = FALSE)
      NULL
    }
  )

  if (is.null(fit)) {
    return(empty)
  }

  coef_names <- names(coef(fit))
  idx <- which(
    str_detect(coef_names, "arm") &
      str_detect(coef_names, "time") &
      str_count(coef_names, ":") == 1
  )

  if (length(idx) != 1) {
    warning("Could not identify arm-by-time coefficient for ", outcome, call. = FALSE)
    return(empty)
  }

  sm <- summary(fit)
  est <- coef(fit)[idx]
  se <- sqrt(diag(vcov(fit)))[idx]
  low <- est - 1.96 * se
  high <- est + 1.96 * se
  p_val <- sm$coefficients[idx, "Pr(>|W|)"]

  list(
    estimate = format_est_ci(est, low, high),
    p_value = p_val
  )
}

get_p_interaction <- function(data, outcome, subgroup_var) {
  data_model <- data %>%
    filter(
      !is.na(.data[[outcome]]),
      !is.na(arm),
      !is.na(time),
      !is.na(.data[[subgroup_var]])
    ) %>%
    droplevels()

  if (
    n_distinct(data_model$arm) < 2 ||
    n_distinct(data_model$time) < 2 ||
    n_distinct(data_model[[subgroup_var]]) < 2
  ) {
    return(NA_real_)
  }

  fit <- tryCatch(
    geeglm(
      as.formula(paste0(outcome, " ~ arm * time * ", subgroup_var)),
      id = id,
      data = data_model,
      family = gaussian(link = "identity"),
      corstr = "exchangeable"
    ),
    error = function(e) {
      warning(
        "Interaction model failed for ", outcome, " by ", subgroup_var,
        ": ", conditionMessage(e),
        call. = FALSE
      )
      NULL
    }
  )

  if (is.null(fit)) {
    return(NA_real_)
  }

  coef_names <- names(coef(fit))
  idx <- which(
    str_detect(coef_names, "arm") &
      str_detect(coef_names, "time") &
      str_detect(coef_names, fixed(subgroup_var)) &
      str_count(coef_names, ":") == 2
  )

  if (length(idx) == 0) {
    warning(
      "Could not identify three-way interaction for ",
      outcome, " by ", subgroup_var,
      call. = FALSE
    )
    return(NA_real_)
  }

  beta <- coef(fit)[idx]
  vcov_matrix <- vcov(fit)[idx, idx, drop = FALSE]
  chisq <- as.numeric(t(beta) %*% safe_solve(vcov_matrix) %*% beta)
  df_test <- length(beta)
  pchisq(chisq, df = df_test, lower.tail = FALSE)
}

args_cmd <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args_cmd, value = TRUE)
script_dir <- if (length(file_arg)) dirname(normalizePath(gsub("~+~", " ", sub("^--file=", "", file_arg[1]), fixed = TRUE))) else {
  getwd()
}
root <- file.path(script_dir, "data-code-table")
project_root <- root

data_path <- file.path(
  project_root,
  "analysis_data",
  "processed",
  "subgroup_analysis_long.xlsx"
)
results_dir <- file.path(script_dir, "analysis_tables")
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(script_dir, "analysis_tables", "details/Table_6"), recursive = TRUE, showWarnings = FALSE)

df <- read_excel(data_path)
df <- df %>%
  mutate(
    id = factor(id),
    arm = factor(arm, levels = c("Control", "Intervention")),
    time = factor(time, levels = c("Baseline", "Followup")),
    region = factor(
      region,
      levels = c("Megacity", "Urban", "Rural")
    ),
    parent_gender = factor(
      parent_gender,
      levels = c("Mother", "Father")
    ),
    age_group = factor(
      ifelse(age_parents >= 40, ">=40", "<40"),
      levels = c("<40", ">=40")
    ),
    higher_education = factor(
      higher_education,
      levels = c("No", "Yes")
    ),
    formal_employment = factor(
      formal_employment,
      levels = c("No", "Yes")
    ),
    household_income = case_when(
      is.na(household_income) ~ NA_character_,
      grepl("^>", trimws(as.character(household_income))) ~ ">200k CNY",
      TRUE ~ "<=200k CNY"
    ),
    household_income = factor(
      household_income,
      levels = c("<=200k CNY", ">200k CNY")
    ),
    hpv_vaccination = factor(
      hpv_vaccination,
      levels = c("No", "Yes")
    ),
    misinformation_exposure = factor(
      misinformation_exposure,
      levels = c("No", "Yes")
    )
  ) %>%
  arrange(id, time)

outcomes <- c(
  module1_score = "Module 1 Coef 95% CI",
  module2_score = "Module 2 Coef 95% CI",
  overall_score = "Overall Coef 95% CI"
)

subgroup_vars <- c(
  "region",
  "parent_gender",
  "age_group",
  "higher_education",
  "formal_employment",
  "household_income",
  "hpv_vaccination",
  "misinformation_exposure"
)

subgroup_labels <- c(
  region = "Region",
  parent_gender = "Parenthood",
  age_group = "Age of parents",
  higher_education = "Higher education",
  formal_employment = "Formal employment",
  household_income = "Household income",
  hpv_vaccination = "HPV vaccination history",
  misinformation_exposure = "Prior misinformation exposure"
)

make_table_structure <- function(data, subgroup_vars, subgroup_labels) {
  overall_row <- tibble(
    row_type = "overall",
    subgroup_var = "Overall",
    level = "",
    Subgroup = "Overall",
    N = n_distinct(data$id)
  )

  subgroup_rows <- map_dfr(subgroup_vars, function(var) {
    label <- subgroup_labels[[var]]
    header_row <- tibble(
      row_type = "subgroup_header",
      subgroup_var = var,
      level = "",
      Subgroup = label,
      N = n_distinct(data$id[!is.na(data[[var]])])
    )

    lvls <- levels(data[[var]])
    lvls <- lvls[!is.na(lvls)]

    level_rows <- map_dfr(lvls, function(lv) {
      dat_sub <- data %>%
        filter(.data[[var]] == lv)

      tibble(
        row_type = "subgroup_level",
        subgroup_var = var,
        level = lv,
        Subgroup = paste0("  ", lv),
        N = n_distinct(dat_sub$id)
      )
    })

    bind_rows(header_row, level_rows)
  })

  bind_rows(overall_row, subgroup_rows)
}

final_table <- make_table_structure(
  data = df,
  subgroup_vars = subgroup_vars,
  subgroup_labels = subgroup_labels
)

p_int_raw <- expand.grid(
  subgroup_var = subgroup_vars,
  outcome = names(outcomes),
  stringsAsFactors = FALSE
) %>%
  as_tibble()

p_int_raw$p_raw <- NA_real_
for (i in seq_len(nrow(p_int_raw))) {
  p_int_raw$p_raw[i] <- get_p_interaction(
    data = df,
    outcome = p_int_raw$outcome[i],
    subgroup_var = p_int_raw$subgroup_var[i]
  )
}

p_int_adjusted <- if (adjust_family == "outcome") {
  p_int_raw %>%
    group_by(outcome) %>%
    mutate(p_adj = p.adjust(p_raw, method = adjust_method)) %>%
    ungroup()
} else {
  p_int_raw %>%
    mutate(p_adj = p.adjust(p_raw, method = adjust_method))
}

lookup_p_adj <- function(outcome, subgroup_var) {
  hit <- p_int_adjusted$p_adj[
    p_int_adjusted$outcome == outcome &
      p_int_adjusted$subgroup_var == subgroup_var
  ]
  if (length(hit) != 1) NA_real_ else hit
}

for (outcome in names(outcomes)) {
  coef_col <- outcomes[[outcome]]
  p_col <- paste0("P int ", str_remove(coef_col, " Coef 95% CI"))

  final_table[[coef_col]] <- ""
  final_table[[p_col]] <- ""

  for (i in seq_len(nrow(final_table))) {
    row_type <- final_table$row_type[i]
    subgroup_var <- final_table$subgroup_var[i]
    level <- final_table$level[i]

    # The overall treatment effect is outside the subgroup correction family.
    if (row_type == "overall") {
      effect <- get_arm_time_effect(
        data = df,
        outcome = outcome
      )

      final_table[[coef_col]][i] <- effect$estimate
      final_table[[p_col]][i] <- format_p(effect$p_value)
    } else if (row_type == "subgroup_header") {
      final_table[[p_col]][i] <- format_p(lookup_p_adj(outcome, subgroup_var))
    } else if (row_type == "subgroup_level") {
      dat_sub <- df %>%
        filter(.data[[subgroup_var]] == level)

      effect <- get_arm_time_effect(
        data = dat_sub,
        outcome = outcome
      )

      final_table[[coef_col]][i] <- effect$estimate
    }
  }
}

final_table_clean <- final_table %>%
  dplyr::select(
    Subgroup,
    N,
    `Module 1 Coef 95% CI`,
    `P int Module 1`,
    `Module 2 Coef 95% CI`,
    `P int Module 2`,
    `Overall Coef 95% CI`,
    `P int Overall`
  )

cat("\nInteraction p-values before and after", adjust_method,
    "correction (family:", adjust_family, ")\n\n")
p_int_adjusted %>%
  mutate(
    outcome = str_remove(outcome, "_score"),
    subgroup = unname(subgroup_labels[subgroup_var]),
    p_raw = formatC(p_raw, format = "f", digits = 4),
    p_adj = formatC(p_adj, format = "f", digits = 4)
  ) %>%
  dplyr::select(outcome, subgroup, p_raw, p_adj) %>%
  as.data.frame() %>%
  print(row.names = FALSE)

p_int_export <- p_int_adjusted %>%
  mutate(
    subgroup = unname(subgroup_labels[subgroup_var]),
    adjust_method = adjust_method,
    adjust_family = adjust_family
  ) %>%
  group_by(outcome) %>%
  mutate(n_tests_in_family = sum(!is.na(p_raw))) %>%
  ungroup() %>%
  dplyr::select(
    outcome, subgroup_var, subgroup, p_raw, p_adj,
    adjust_method, adjust_family, n_tests_in_family
  )

if (adjust_family != "outcome") {
  p_int_export$n_tests_in_family <- sum(!is.na(p_int_export$p_raw))
}

out_p_csv <- file.path(results_dir, "details/Table_6/subgroup_interaction_p_multiplicity.csv")
write.csv(p_int_export, out_p_csv, row.names = FALSE)
cat("\nSaved:", out_p_csv, "\n")

out_csv <- file.path(
  results_dir,
  "Table6.csv"
)
write.csv(final_table_clean, out_csv, row.names = FALSE)
cat("\nSaved:", out_csv, "\n")
