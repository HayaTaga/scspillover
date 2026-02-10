# ==============================================================
# State tobacco panel + state adjacency (rook) with CA vector
# Sources:
#   - Smoking panel: https://github.com/scunning1975/mixtape
#   - State boundaries: U.S. Census Bureau, TIGER/Line (tl_2024_us_state)
# ==============================================================

# --- Libraries
library(haven) # read_dta
library(dplyr)
library(readr)
library(sf)
library(spdep)
library(tibble) # rownames_to_column

# --------------------------------------------------------------
# 1) Load & tidy the smoking panel (state, state_id, year, cigsale, retprice)
# --------------------------------------------------------------
panel <- read_dta("data/raw/smoking.dta") %>%
  mutate(
    # Preserve Stata numeric code as state_id
    state_id = as.integer(as.numeric(state)),
    # Convert labelled state to human-readable string
    state = as.character(as_factor(state))
  ) %>%
  select(state, state_id, year, cigsale, retprice)

# Persist the tidy panel for reproducibility
write_csv(panel, "data/processed/tobacco/us_state_panel.csv")

# --------------------------------------------------------------
# 2) Load state polygons and create rook-contiguity adjacency
# --------------------------------------------------------------
# Read the state shapefile (states + DC; territories might be present depending on source)
sf_obj <- st_read(
  "data/raw/tl_2024_us_state/tl_2024_us_state.shp",
  quiet = TRUE
)

# Keep only the 50 states + DC (drop Puerto Rico & territories)
# STUSPS is the two-letter USPS code
territories_drop <- c("PR", "VI", "GU", "MP", "AS")
sf_obj <- sf_obj %>%
  filter(!(STUSPS %in% territories_drop))

# Ensure valid geometries (robustness for topology operations)
sf_obj <- sf_obj %>% mutate(geometry = st_make_valid(geometry))

# Build neighbor list using rook contiguity (shared border; no corner-touch)
# poly2nb works directly on sf objects
nb_all <- poly2nb(sf_obj, queen = FALSE)

# Binary adjacency matrix (0/1). Use NAME as dimension names
W_all <- nb2mat(nb_all, style = "B", zero.policy = TRUE)
rownames(W_all) <- colnames(W_all) <- sf_obj$NAME

# --------------------------------------------------------------
# 3) Restrict adjacency to the states present in the panel
#    (Intersection; keep panel names order for downstream consistency)
# --------------------------------------------------------------
states_panel <- panel %>% distinct(state) %>% pull() %>% as.character()

# Intersect by name; warn on mismatches
missing_in_map <- setdiff(states_panel, rownames(W_all))
missing_in_panel <- setdiff(rownames(W_all), states_panel)
if (length(missing_in_map) > 0) {
  warning(
    "States in panel but not in map (will be dropped): ",
    paste(missing_in_map, collapse = ", ")
  )
}
# Subset adjacency by the states that exist in both
states_keep <- intersect(states_panel, rownames(W_all))

# Re-index adjacency to the chosen set, preserving the panel order
idx_keep <- match(states_keep, rownames(W_all))
W_sub <- W_all[idx_keep, idx_keep, drop = FALSE]
rownames(W_sub) <- colnames(W_sub) <- states_keep

# Optional sanity: symmetry and diagonal
stopifnot(all(W_sub == t(W_sub)))
stopifnot(all(diag(W_sub) == 0))

# --------------------------------------------------------------
# 4) Extract California’s adjacency vector, then drop California
# --------------------------------------------------------------
ca_name <- "California"
if (!ca_name %in% rownames(W_sub)) {
  stop("California is not in the restricted state set.")
}
# Vector of CA adjacency against the restricted set (including CA itself = 0)
ca_vec <- W_sub[ca_name, , drop = TRUE]

# Drop California from the matrix (row & column)
keep_no_ca <- setdiff(rownames(W_sub), ca_name)
W_no_ca <- W_sub[keep_no_ca, keep_no_ca, drop = FALSE]

# Drop California from the vector and align order with W_no_ca columns
ca_vec_no_ca <- ca_vec[keep_no_ca]

# --------------------------------------------------------------
# 5) Persist outputs (friendly, rectangular CSVs)
# --------------------------------------------------------------
# (a) Full rook adjacency for the restricted set
W_sub_df <- W_sub %>%
  as.data.frame() %>%
  rownames_to_column(var = "state")
write_csv(W_sub_df, "data/processed/tobacco/adjacency_rook_panel.csv")

# (b) Rook adjacency with California removed
W_no_ca_df <- W_no_ca %>%
  as.data.frame() %>%
  rownames_to_column(var = "state")
write_csv(W_no_ca_df, "data/processed/tobacco/weight_mat.csv")

# (c) California adjacency vector (matching W_no_ca column order)
ca_vec_df <- tibble(
  state = names(ca_vec_no_ca),
  adj = as.numeric(ca_vec_no_ca)
)
write_csv(ca_vec_df, "data/processed/tobacco/weight_vec.csv")

# --------------------------------------------------------------
# 6) Minimal diagnostics (printed)
# --------------------------------------------------------------
message(sprintf(
  "Panel states: %d, Map states kept: %d, Matrix dim: %dx%d",
  length(states_panel),
  length(states_keep),
  nrow(W_sub),
  ncol(W_sub)
))
message(sprintf(
  "Matrix without CA: %dx%d; CA vector length: %d",
  nrow(W_no_ca),
  ncol(W_no_ca),
  length(ca_vec_no_ca)
))


rm(list = ls())
panel_df <- read_csv("data/processed/tobacco/us_state_panel.csv")
w <- read_csv("data/processed/tobacco/weight_vec.csv")
W <- read_csv("data/processed/tobacco/weight_mat.csv")

california_smoking <- list(
  panel_df = panel_df,
  w = w,
  W = W
)

usethis::use_data(california_smoking, overwrite = TRUE)
