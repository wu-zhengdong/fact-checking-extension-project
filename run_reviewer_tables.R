# Reproduce Table 2 and Supplementary Tables S2, S3, S4, S7 and S8.
# Usage: Rscript run_reviewer_tables.R [path to deidentified_item_level_responses.csv]

a <- grep("^--file=", commandArgs(FALSE), value = TRUE)
here <- if (length(a)) dirname(normalizePath(gsub("~+~", " ", sub("^--file=", "", a[1]), fixed = TRUE))) else {
  getwd()
}
source(file.path(here, "lib", "reviewer_data.R"))
root <- here
data_path <- reviewer_data_path(root)
message("Data file: ", data_path)

scripts <- c("Table 2.R", paste0("Supplementary Table ", c("S2", "S3", "S4", "S7", "S8"), ".R"))
for (script in scripts) {
  message("Running ", script)
  status <- system2(file.path(R.home("bin"), "Rscript"),
                    shQuote(c(file.path(here, script), data_path)))
  if (status != 0L) stop(script, " failed with exit status ", status)
}

out <- file.path(root, "analysis_tables")
dir.create(file.path(out, "details"), recursive = TRUE, showWarnings = FALSE)
# The tables ran in separate processes, so attach the packages for sessionInfo().
pkgs <- c("dplyr", "tidyr", "geepack", "lme4")
for (p in pkgs) suppressPackageStartupMessages(library(p, character.only = TRUE))
writeLines(c(capture.output(sessionInfo()), "", "Analysis package versions:",
             sprintf("  %s %s", pkgs,
                     vapply(pkgs, function(p) as.character(packageVersion(p)), character(1)))),
           file.path(out, "details", "reviewer_sessionInfo.txt"))
message("All six tables saved in ", out)
