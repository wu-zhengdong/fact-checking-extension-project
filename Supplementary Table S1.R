# Supplementary Table S1. Experimental posts, labels and assessment statements.

a <- grep("^--file=", commandArgs(FALSE), value = TRUE)
here <- if (length(a)) dirname(normalizePath(gsub("~+~", " ", sub("^--file=", "", a[1]), fixed = TRUE))) else {
  getwd()
}
x <- data.frame(
`Item` = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11"),
`Module` = c("Module 1", "Module 1", "Module 1", "Module 1", "Module 1", "Module 1", "Module 2", "Module 2", "Module 2", "Module 2", "Module 3"),
`Post topic` = c("Infertility", "Vaccination age", "Vaccine ingredients", "HPV infection", "Antibiotic use", "High-risk HPV types", "Pregnancy/breastfeeding", "Adverse reactions", "Cervical cancer screening", "Vaccination timing", "Reproductive health"),
`Post veracity` = c("Misinformation", "Accurate", "Misinformation", "Accurate", "Misinformation", "Accurate", "Accurate", "Misinformation", "Accurate", "Misinformation", "Misinformation"),
`Fact-check label` = c("Red", "Green", "Red", "Green", "Red", "Green", "N/A", "N/A", "N/A", "N/A", "Red if selected"),
`Information-assessment statement` = c("HPV vaccination may cause infertility.", "Girls aged 9–15 years should be prioritized for HPV vaccination.", "HPV vaccination may cause poisoning.", "HPV infection is equivalent to having a “cold” of the cervix.", "HPV vaccination should not be administered while taking antibiotics.", "The most common high-risk HPV types are HPV 16 and 18.", "Women who are breastfeeding should not receive HPV vaccination because vaccine components may enter breast milk and harm the infant.", "HPV vaccination frequently causes severe allergic reactions.", "Cervical cancer screening is no longer necessary after HPV vaccination.", "Bivalent HPV vaccine is enough, early vaccination is more important.", "HPV vaccination may cause premature ovarian failure."),
`Correct answer` = c("FALSE", "TRUE", "FALSE", "FALSE", "FALSE", "TRUE", "FALSE", "FALSE", "FALSE", "TRUE", "FALSE"),
`Label–item alignment` = c("Concordant", "Concordant", "Concordant", "Discordant", "Concordant", "Concordant", "N/A", "N/A", "N/A", "N/A", "Concordant")
, check.names = FALSE)
out <- file.path(here, "analysis_tables")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
write.csv(x, file.path(out, "Supplementary_Table_S1.csv"), row.names = FALSE)
