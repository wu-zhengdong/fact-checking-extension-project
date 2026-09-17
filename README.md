# fact-checking-extension-project

Analysis code for the manuscript *A randomized experiment of a simulated AI-attributed
fact-checking interface for HPV vaccine information assessment in China*. One R script per table,
named after the table it produces.

## Reproducing Table 2 and Supplementary Tables S2, S3, S4, S7 and S8

These six run from `deidentified_item_level_responses.csv`, supplied with the manuscript through the
journal submission system. Download this code with **Code ▸ Download ZIP** or `git clone`, then
install R (tested on 4.6.1) and four packages:

```r
install.packages(c("dplyr", "tidyr", "geepack", "lme4"))
```

From this folder, give the path to the CSV:

```bash
Rscript run_reviewer_tables.R /path/to/deidentified_item_level_responses.csv
```

About 15 seconds. Each table is printed to the console and written to `analysis_tables/`. Compare
them with the corresponding tables in the manuscript.

To run one table at a time, quote the filename, because these filenames contain spaces:

```bash
Rscript "Table 2.R" /path/to/deidentified_item_level_responses.csv
Rscript "Supplementary Table S1.R"     # fixed design table, needs no data
```

From R or RStudio instead of a terminal:

```r
setwd("/path/to/this/folder")
Sys.setenv(FC_DATA = "/path/to/deidentified_item_level_responses.csv")
source("run_reviewer_tables.R")
```

`lib/reviewer_data.R` is the shared reader the six scripts `source()`; keep it where it is.

## License

MIT, see [LICENSE](LICENSE).
