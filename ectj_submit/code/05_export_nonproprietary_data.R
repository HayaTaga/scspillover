dir.create("data/nonproprietary", recursive = TRUE, showWarnings = FALSE)

extract_panel <- function(obj, label) {
  if ("panel" %in% names(obj)) {
    return(obj[["panel"]])
  }
  if ("panel_df" %in% names(obj)) {
    return(obj[["panel_df"]])
  }
  stop(sprintf("Missing panel/panel_df in %s", label))
}

load("data/california_smoking.rda")
write.csv(
  extract_panel(california_smoking, "california_smoking"),
  "data/nonproprietary/california_panel.csv",
  row.names = FALSE
)
write.csv(
  california_smoking$w,
  "data/nonproprietary/california_w_vector.csv",
  row.names = FALSE
)
write.csv(
  california_smoking$W,
  "data/nonproprietary/california_W_matrix.csv",
  row.names = FALSE
)

load("data/sudan_secession.rda")
write.csv(
  extract_panel(sudan_secession, "sudan_secession"),
  "data/nonproprietary/sudan_panel.csv",
  row.names = FALSE
)
write.csv(
  sudan_secession$w,
  "data/nonproprietary/sudan_w_vector.csv",
  row.names = FALSE
)
write.csv(
  sudan_secession$W,
  "data/nonproprietary/sudan_W_matrix.csv",
  row.names = FALSE
)

message("Exported CSV copies into data/nonproprietary/.")
