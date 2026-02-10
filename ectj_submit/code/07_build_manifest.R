classify_file <- function(path) {
  if (grepl("^(code/|R/|src/|vignettes/)", path) || path %in% c("Makefile")) {
    return("code")
  }
  if (grepl("^data/", path) && grepl("[.](rda|csv)$", path)) {
    return("input_data")
  }
  if (grepl("^output/", path)) {
    return("generated_output")
  }
  return("documentation")
}

bundle_flag <- function(path) {
  if (grepl("^output/figures/", path) && basename(path) != ".gitkeep") return("no")
  if (grepl("^output/tables/", path) && basename(path) != ".gitkeep") return("no")
  return("yes")
}

gen_flag <- function(path) {
  if (grepl("^output/figures/", path) && basename(path) != ".gitkeep") return("yes")
  if (grepl("^output/tables/", path) && basename(path) != ".gitkeep") return("yes")
  return("no")
}

files <- list.files(".", recursive = TRUE, all.files = TRUE, no.. = TRUE)
files <- files[!grepl("^\\.git/", files)]
files <- files[file.info(files)$isdir %in% c(FALSE)]
files <- sort(files)

manifest_actual <- data.frame(
  path = files,
  class = vapply(files, classify_file, character(1)),
  present_in_repo = "yes",
  bundled_at_submission = vapply(files, bundle_flag, character(1)),
  generated_by_editor = vapply(files, gen_flag, character(1)),
  source = "actual_file",
  stringsAsFactors = FALSE
)

if (file.exists("RESULTS_MAPPING.csv")) {
  rm_map <- read.csv("RESULTS_MAPPING.csv", stringsAsFactors = FALSE)
  output_specs <- unique(rm_map$output_file)

  exists_for_spec <- function(spec) {
    matches <- Sys.glob(spec)
    if (length(matches) > 0) "yes" else "no"
  }

  manifest_spec <- data.frame(
    path = output_specs,
    class = "generated_output",
    present_in_repo = vapply(output_specs, exists_for_spec, character(1)),
    bundled_at_submission = "no",
    generated_by_editor = "yes",
    source = "results_mapping_spec",
    stringsAsFactors = FALSE
  )

  manifest <- rbind(manifest_actual, manifest_spec)
} else {
  manifest <- manifest_actual
}

manifest <- manifest[order(manifest$path, manifest$source), ]
write.csv(manifest, "MANIFEST.csv", row.names = FALSE)
cat(sprintf("Wrote MANIFEST.csv with %d rows\n", nrow(manifest)))
