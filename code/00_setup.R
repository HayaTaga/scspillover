required_pkgs <- c(
  "Rcpp", "RcppArmadillo", "Matrix", "progress",
  "dplyr", "tidyr", "ggplot2", "stringr",
  "knitr", "rmarkdown", "kableExtra", "ggpattern",
  "purrr", "readr", "glue", "coda", "mcmcse", "posterior"
)

missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))]
if (length(missing_pkgs) > 0) {
  install_missing <- identical(Sys.getenv("SCSPILL_INSTALL_MISSING", "false"), "true")
  if (!install_missing) {
    stop(
      sprintf(
        paste(
          "Missing required packages:",
          "%s",
          "Set SCSPILL_INSTALL_MISSING=true to auto-install from CRAN.",
          sep = "\n"
        ),
        paste(missing_pkgs, collapse = ", ")
      )
    )
  }
  install.packages(missing_pkgs, repos = "https://cloud.r-project.org")
}

# Re-check after optional installation.
still_missing <- required_pkgs[!vapply(required_pkgs, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))]
if (length(still_missing) > 0) {
  stop(sprintf("Could not load required packages: %s", paste(still_missing, collapse = ", ")))
}

lock_path <- "DEPENDENCY_LOCK.csv"
if (file.exists(lock_path)) {
  lock_df <- utils::read.csv(lock_path, stringsAsFactors = FALSE)
  lock_pkg <- subset(lock_df, type == "package")
  lock_pkg <- lock_pkg[lock_pkg$name %in% required_pkgs, , drop = FALSE]
  if (nrow(lock_pkg) > 0) {
    mismatches <- character(0)
    for (i in seq_len(nrow(lock_pkg))) {
      pkg <- lock_pkg$name[i]
      locked_ver <- lock_pkg$version[i]
      installed_ver <- as.character(utils::packageVersion(pkg))
      if (!identical(installed_ver, locked_ver)) {
        mismatches <- c(
          mismatches,
          sprintf("%s (installed=%s, lock=%s)", pkg, installed_ver, locked_ver)
        )
      }
    }
    if (length(mismatches) > 0) {
      enforce_lock <- identical(Sys.getenv("SCSPILL_ENFORCE_LOCK", "false"), "true")
      msg <- paste(
        "Version differences against DEPENDENCY_LOCK.csv:",
        paste(mismatches, collapse = "; "),
        sep = "\n"
      )
      if (enforce_lock) {
        stop(msg)
      } else {
        message(msg)
      }
    }
  }
}

pkg_versions <- vapply(required_pkgs, function(p) as.character(utils::packageVersion(p)), character(1))
message("All required packages are available.")
message(sprintf("R version: %s", R.version.string))
message(sprintf("Platform: %s", paste(Sys.info()[c("sysname", "release", "machine")], collapse = " ")))
message("Package versions:")
for (p in names(pkg_versions)) {
  message(sprintf("- %s: %s", p, pkg_versions[[p]]))
}
