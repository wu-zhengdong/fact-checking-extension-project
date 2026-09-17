# Table 4. Voluntary adoption of the extension.
# Log-binomial risk ratio comparing randomized arms.

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
dir.create(file.path(script_dir, "analysis_tables", "details/Table_4"), recursive = TRUE, showWarnings = FALSE)

OUTCOME_LABEL <- "Chose to use the extension"

fmt_n_pct <- function(n, denom) {
  sprintf("%s (%.1f)", format(n, big.mark = ","), 100 * n / denom)
}

fmt_effect <- function(rr, lo, hi) {
  sprintf("%.2f (%.2f-%.2f)", rr, lo, hi)
}

fmt_p <- function(p) {
  if (is.na(p)) return("")
  if (p < 0.001) return("<0.001")
  sprintf("%.3f", p)
}

# Adoption is measured once; use one post-exposure row per participant.
panel <- read_xlsx(data_path, sheet = "panel")

varying <- panel %>%
  group_by(pid) %>%
  summarise(k = n_distinct(switch_3_1), .groups = "drop")
if (any(varying$k != 1)) stop("switch_3_1 varies within a participant")

dat <- panel %>%
  filter(time == 1) %>%
  mutate(
    adopt = as.integer(switch_3_1 == 1),
    class_int = as.integer(factor(class_id))
  ) %>%
  arrange(class_int, pid)

if (any(is.na(dat$class_id))) stop("Missing class_id values")
stopifnot(
  nrow(dat) == 1548L,
  !anyNA(dat$adopt),
  sum(dat$treatment == 1) == 771L,
  sum(dat$treatment == 0) == 777L
)

one_class_per_pid <- dat %>%
  distinct(pid, class_id) %>%
  count(pid, name = "n")
if (any(one_class_per_pid$n != 1)) stop("A participant maps to >1 class_id")

n_total <- nrow(dat)
n_fc <- sum(dat$treatment == 1)
n_sv <- sum(dat$treatment == 0)
n_classes <- n_distinct(dat$class_int)

k_total <- sum(dat$adopt)
k_fc <- sum(dat$adopt == 1 & dat$treatment == 1)
k_sv <- sum(dat$adopt == 1 & dat$treatment == 0)
stopifnot(k_total == 1100L, k_fc == 563L, k_sv == 537L, n_classes == 218L)

cells <- list(
  total = fmt_n_pct(k_total, n_total),
  fc = fmt_n_pct(k_fc, n_fc),
  sv = fmt_n_pct(k_sv, n_sv)
)

fit_glm <- glm(
  adopt ~ treatment,
  family = binomial(link = "log"),
  data = dat
)
if (!fit_glm$converged) stop("Log-binomial GLM did not converge")

fit_gee <- geeglm(
  adopt ~ treatment,
  id = class_int,
  data = dat,
  family = binomial(link = "log"),
  corstr = "independence"
)

extract <- function(estimate, se, p) {
  list(
    log_rr = estimate,
    se = se,
    rr = exp(estimate),
    lo = exp(estimate - 1.96 * se),
    hi = exp(estimate + 1.96 * se),
    p = p
  )
}

co_glm <- summary(fit_glm)$coefficients["treatment", ]
co_gee <- summary(fit_gee)$coefficients["treatment", ]

est_glm <- extract(co_glm[["Estimate"]], co_glm[["Std. Error"]], co_glm[["Pr(>|z|)"]])
est_gee <- extract(co_gee[["Estimate"]], co_gee[["Std.err"]], co_gee[["Pr(>|W|)"]])

stopifnot(abs(est_glm$log_rr - est_gee$log_rr) < 1e-8)
stopifnot(abs(est_glm$rr - (k_fc / n_fc) / (k_sv / n_sv)) < 1e-8)

specs <- list(
  list(
    key = "main",
    model = "Log-binomial GLM (participants independent)",
    se_type = "Model-based Wald standard error",
    clusters = NA_integer_,
    est = est_glm
  ),
  list(
    key = "class_cluster",
    model = "Log-binomial GEE, class-cluster-robust standard errors",
    se_type = "Cluster-robust (sandwich) standard error, school class",
    clusters = n_classes,
    est = est_gee
  )
)

table_rows <- lapply(specs, function(spec) {
  data.frame(
    Model = spec$model,
    `Standard error` = spec$se_type,
    Outcome = OUTCOME_LABEL,
    Total = cells$total,
    `Fact-checking group` = cells$fc,
    `Standard-view group` = cells$sv,
    `Effect estimate (95% CI)` = fmt_effect(
      spec$est$rr, spec$est$lo, spec$est$hi
    ),
    `P value` = fmt_p(spec$est$p),
    `Total N` = n_total,
    `Fact-checking group N` = n_fc,
    `Standard-view group N` = n_sv,
    `Number of classes` = spec$clusters,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
})

out_table <- bind_rows(table_rows)
out_csv <- file.path(out_dir, "Table4.csv")
write.csv(out_table[1, ], out_csv, row.names = FALSE)

chisq_p <- suppressWarnings(
  chisq.test(table(dat$treatment, dat$adopt), correct = FALSE)$p.value
)

use_lme4 <- requireNamespace("lme4", quietly = TRUE)
if (use_lme4) {
  icc_fit <- lme4::lmer(adopt ~ 1 + (1 | class_id), data = dat, REML = TRUE)
  vc <- as.data.frame(lme4::VarCorr(icc_fit))
  var_class <- vc$vcov[vc$grp == "class_id"]
  var_resid <- vc$vcov[vc$grp == "Residual"]
} else {
  icc_fit <- nlme::lme(
    adopt ~ 1, random = ~ 1 | class_id, data = dat, method = "REML"
  )
  sds <- as.numeric(nlme::VarCorr(icc_fit)[, "StdDev"])
  var_class <- sds[1]^2
  var_resid <- sds[2]^2
}
icc <- var_class / (var_class + var_resid)

class_sizes <- dat %>% count(class_int, name = "n")

details <- bind_rows(lapply(specs, function(spec) {
  est <- spec$est
  se_corrected <- if (is.na(spec$clusters)) {
    NA_real_
  } else {
    est$se * sqrt(spec$clusters / (spec$clusters - 1))
  }
  data.frame(
    Model = spec$model,
    `Standard error` = spec$se_type,
    Formula = "adopt ~ treatment",
    Family = "binomial(link = \"log\")",
    Outcome = OUTCOME_LABEL,
    `log RR` = sprintf("%.6f", est$log_rr),
    `SE (log RR)` = sprintf("%.6f", est$se),
    `SE with G/(G-1) correction` = if (is.na(se_corrected)) {
      ""
    } else {
      sprintf("%.6f", se_corrected)
    },
    `SE ratio vs main` = sprintf("%.4f", est$se / est_glm$se),
    RR = sprintf("%.4f", est$rr),
    `95% CI lower` = sprintf("%.4f", est$lo),
    `95% CI upper` = sprintf("%.4f", est$hi),
    `P value` = sprintf("%.4f", est$p),
    `Adoption risk, fact-checking` = sprintf("%.4f", k_fc / n_fc),
    `Adoption risk, standard view` = sprintf("%.4f", k_sv / n_sv),
    `Risk difference` = sprintf("%.4f", k_fc / n_fc - k_sv / n_sv),
    `Pearson chi-square P (reference)` = sprintf("%.4f", chisq_p),
    `Class-level ICC of adoption` = sprintf("%.4f", icc),
    `Number of classes` = n_classes,
    `Participants per class, median` = stats::median(class_sizes$n),
    `Participants per class, range` = sprintf(
      "%d-%d", min(class_sizes$n), max(class_sizes$n)
    ),
    N = n_total,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}))

out_details <- file.path(out_dir, "details/Table_4/table4_adoption_risk_ratio_details.csv")
write.csv(details, out_details, row.names = FALSE)

cat(strrep("=", 70), "\n", sep = "")
cat("Behavioral adoption of the extension as an outcome\n")
cat(strrep("=", 70), "\n", sep = "")
cat(sprintf(
  "N = %s (%s fact-checking, %s standard view) in %d school classes\n",
  format(n_total, big.mark = ","), format(n_fc, big.mark = ","),
  format(n_sv, big.mark = ","), n_classes
))
cat(sprintf(
  "Adopted: %s total, %s fact-checking, %s standard view\n\n",
  cells$total, cells$fc, cells$sv
))

for (spec in specs) {
  est <- spec$est
  cat(spec$model, "\n")
  cat(sprintf("  %s\n", spec$se_type))
  cat(sprintf(
    "  RR %.4f (95%% CI %.4f-%.4f), SE(log RR) %.5f, P %.4f\n\n",
    est$rr, est$lo, est$hi, est$se, est$p
  ))
}

cat(sprintf(
  "SE(log RR) %.5f -> %.5f (x%.3f) when clustering on class\n",
  est_glm$se, est_gee$se, est_gee$se / est_glm$se
))
cat(sprintf("Class-level ICC of adoption: %.4f (%s)\n",
            icc, if (use_lme4) "lme4::lmer" else "nlme::lme"))
cat(sprintf("Pearson chi-square P for reference (table4.R): %s\n\n",
            fmt_p(chisq_p)))

print(
  out_table[, c(
    "Model", "Outcome", "Total", "Fact-checking group",
    "Standard-view group", "Effect estimate (95% CI)", "P value"
  )],
  row.names = FALSE
)

cat("\nSaved:", out_csv, "\n")
cat("Saved:", out_details, "\n")
