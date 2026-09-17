# Run all 17 tables. Tables 1, 3, 4, 5, 6 and Supplementary Tables S5, S6, S9, S10, S11
# additionally require the restricted workbooks described in the README, so this runner
# only completes for users who hold those files.

a <- grep("^--file=", commandArgs(FALSE), value = TRUE)
here <- if (length(a)) dirname(normalizePath(gsub("~+~", " ", sub("^--file=", "", a[1]), fixed = TRUE))) else {
  getwd()
}
extra <- commandArgs(trailingOnly = TRUE)
scripts <- c(paste0("Table ", 1:6, ".R"), paste0("Supplementary Table S", 1:11, ".R"))
for (script in scripts) {
  message("Running ", script)
  status <- system2(file.path(R.home("bin"), "Rscript"),
                    shQuote(c(file.path(here, script), extra)))
  if (status != 0L) stop(script, " failed with exit status ", status)
}

out <- file.path(here, "analysis_tables")
dir.create(file.path(out, "details"), recursive = TRUE, showWarnings = FALSE)
# Each table ran in its own Rscript process, so attach the packages here to make
# sessionInfo() report the versions that actually produced the tables.
pkgs <- c("dplyr", "tidyr", "geepack", "lme4", "readxl", "purrr", "stringr", "tibble", "MASS")
pkgs <- pkgs[vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
for (p in pkgs) suppressPackageStartupMessages(library(p, character.only = TRUE))
writeLines(c(capture.output(sessionInfo()), "", "Analysis package versions:",
             sprintf("  %s %s", pkgs,
                     vapply(pkgs, function(p) as.character(packageVersion(p)), character(1)))),
           file.path(out, "details", "internal_sessionInfo.txt"))
