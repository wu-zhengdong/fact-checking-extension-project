# Read and validate the de-identified item-response data.

# The CSV is distributed through the journal submission system, not with this code, so
# accept it from a command-line argument, the FC_DATA variable, or this folder.
reviewer_data_path <- function(repo_root) {
  name <- "deidentified_item_level_responses.csv"
  candidates <- c(commandArgs(trailingOnly = TRUE), Sys.getenv("FC_DATA"),
                  file.path(repo_root, name), file.path(repo_root, "data", name), name)
  candidates <- candidates[nzchar(candidates)]
  for (p in candidates) {
    if (dir.exists(p)) p <- file.path(p, name)
    if (file.exists(p) && !dir.exists(p)) return(normalizePath(p))
  }
  stop("Cannot find ", name, ".\nGive its path, for example\n",
       "  Rscript run_reviewer_tables.R /path/to/", name,
       "\nor set FC_DATA, or copy the file into ", repo_root,
       "\nLooked in:\n", paste0("  ", candidates, collapse = "\n"), call. = FALSE)
}

read_reviewer_items <- function(path) {
  # Preserve TRUE/FALSE response labels as text.
  x <- read.csv(path, colClasses = "character", check.names = FALSE)
  expected <- c("participant_id", "group", "time", "item_id", "analysis_module",
                "post_veracity", "correct_answer", "response", "accuracy")
  stopifnot(identical(names(x), expected), !anyNA(x))
  stopifnot(all(x$item_id %in% as.character(1:10)), all(x$accuracy %in% c("0", "1")))
  x$item_id <- as.integer(x$item_id)
  x$accuracy <- as.integer(x$accuracy)
  stopifnot(nrow(x) == 30960L, length(unique(x$participant_id)) == 1548L,
            all(grepl("^P[0-9]{4}$", x$participant_id)),
            all(x$group %in% c("fact-checking", "standard-view")),
            all(x$time %in% c("pre", "post")),
            all(x$response %in% c("TRUE", "FALSE", "Don't know")),
            !anyDuplicated(x[c("participant_id", "time", "item_id")]),
            all(table(x$participant_id, x$time) == 10L),
            all(x$accuracy == as.integer(x$response == x$correct_answer)))
  stopifnot(all(vapply(split(x$group, x$participant_id), function(z) length(unique(z)), integer(1)) == 1L))
  arms <- table(unique(x[c("participant_id", "group")])$group)
  stopifnot(arms[["fact-checking"]] == 771L, arms[["standard-view"]] == 777L)
  key <- c("FALSE", "TRUE", "FALSE", "FALSE", "FALSE", "TRUE", "FALSE", "FALSE", "FALSE", "TRUE")
  veracity <- c("misinformation", "accurate", "misinformation", "accurate", "misinformation", "accurate", "accurate", "misinformation", "accurate", "misinformation")
  stopifnot(all(x$correct_answer == key[x$item_id]),
            all(x$post_veracity == veracity[x$item_id]),
            all(x$analysis_module == ifelse(x$item_id <= 6L, "Module 1", "Module 2")))
  x[order(x$participant_id, match(x$time, c("pre", "post")), x$item_id), ]
}

reviewer_panel <- function(path) {
  x <- read_reviewer_items(path)
  # S1 item numbers mapped to statements* variables; pid remains the anonymous ID.
  source_index <- c(8L, 2L, 9L, 4L, 5L, 1L, 6L, 10L, 7L, 3L)
  panel <- unique(x[c("participant_id", "group", "time")])
  names(panel)[1] <- "pid"
  panel$treatment <- as.integer(panel$group == "fact-checking")
  panel$time <- as.integer(panel$time == "post")
  panel$group <- NULL
  for (item in 1:10) {
    z <- x[x$item_id == item, ]
    stopifnot(identical(z$participant_id, panel$pid),
              identical(as.integer(z$time == "post"), panel$time))
    var <- paste0("statements", source_index[item])
    panel[[var]] <- match(z$response, c("TRUE", "FALSE", "Don't know"))
    panel[[paste0(var, "_s")]] <- z$accuracy
  }
  panel
}
