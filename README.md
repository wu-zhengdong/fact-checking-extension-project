# fact-checking-extension-project

Analysis code for the manuscript *A randomized experiment of a simulated AI-attributed
fact-checking interface for HPV vaccine information*.

## Reproducing Table 2 and Supplementary Tables S2, S3, S4, S7 and S8

These six tables use `deidentified_item_level_responses.csv`, provided confidentially through the
journal submission system. With R installed (tested on 4.6.1), install the required packages:

```r
install.packages(c("dplyr", "tidyr", "geepack", "lme4"))
```

From the repository folder, run:

```bash
Rscript run_reviewer_tables.R /path/to/deidentified_item_level_responses.csv
```

Results are saved in `analysis_tables/`. Table 2 includes Cohen's d and its 95% confidence interval.
Other statistical tables require additional study data not included in this CSV.

To run Table 2 alone:

```bash
Rscript "Table 2.R" /path/to/deidentified_item_level_responses.csv
```

From R or RStudio instead of a terminal:

```r
setwd("/path/to/this/folder")
Sys.setenv(FC_DATA = "/path/to/deidentified_item_level_responses.csv")
source("run_reviewer_tables.R")
```

## License

MIT, see [LICENSE](LICENSE).
