# Evaluate NCT split conditions without executing arbitrary R code.
evaluate_nct_filter <- function(filter, data) {
  if (length(filter) != 1L || is.na(filter) ||
      trimws(filter) %in% c("", "NA")) {
    return(rep(TRUE, nrow(data)))
  }

  conditions <- strsplit(filter, "&", fixed = TRUE)[[1]]
  keep <- rep(TRUE, nrow(data))

  pattern <- paste0(
    "^(?:data_tree\\$|data_est\\$)?",
    "(X\\.[0-9]+)\\s*",
    "(>=|<=|==|>|<)\\s*",
    "(-?(?:[0-9]+\\.?[0-9]*|\\.[0-9]+)(?:[eE][+-]?[0-9]+)?)$"
  )

  for (condition in conditions) {
    condition <- trimws(condition)
    match <- regmatches(
      condition,
      regexec(pattern, condition, perl = TRUE)
    )[[1]]

    if (length(match) != 4L) {
      stop("Unsupported NCT filter condition: ", condition)
    }

    variable <- match[2]
    operator <- match[3]
    threshold <- as.numeric(match[4])

    if (!variable %in% names(data)) {
      stop("Unknown NCT filter variable: ", variable)
    }

    x <- data[[variable]]

    comparison <- switch(
      operator,
      ">=" = x >= threshold,
      "<=" = x <= threshold,
      "==" = x == threshold,
      ">" = x > threshold,
      "<" = x < threshold
    )

    comparison[is.na(comparison)] <- FALSE
    keep <- keep & comparison
  }

  keep
}
