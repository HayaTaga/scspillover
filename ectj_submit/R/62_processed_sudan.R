# ============================================
# Replication script: Africa panel (2000–2015)
# Preprocessing with Sudan/South Sudan rules
# ============================================

# --- Libraries
library(tidyverse) # dplyr, tidyr, stringr, purrr, readr, etc.
library(httr)
library(jsonlite)

# -------------------------------
# 0) Load raw WDI-like CSV
# -------------------------------
raw <- read.csv(
  "data/raw/P_Data_Extract_From_World_Development_Indicators/24cc87ae-8d91-47e6-b829-e467008e1958_Data.csv",
  na = c("..", ""),
  stringsAsFactors = FALSE
)

# Detect year columns safely: patterns like "X2000..YR2000." or "2000 [YR2000]"
nm <- names(raw)
year_cols <- nm[stringr::str_detect(nm, "^X?(19|20)\\d{2}")]
if (length(year_cols) == 0) {
  stop("No year columns detected.")
}

# Long format (country x series x year)
long_dat <- raw %>%
  tidyr::pivot_longer(
    cols = dplyr::all_of(year_cols),
    names_to = "year_raw",
    values_to = "value"
  ) %>%
  dplyr::mutate(
    year = as.integer(stringr::str_match(year_raw, "^X?(\\d{4})")[, 2])
  ) %>%
  dplyr::select(
    Country.Name,
    Country.Code,
    Series.Name,
    Series.Code,
    year,
    value
  )

# Keep only required series (as per your mapping)
target_codes <- c(
  "NY.GDP.PCAP.KD", # GDP per capita (constant 2015 US$)
  "NY.GDP.MKTP.KD", # GDP (constant 2015 US$)
  "SP.POP.TOTL", # Population, total
  "NE.EXP.GNFS.ZS", # Exports of goods and services (% of GDP)
  "TG.VAL.TOTL.GD.ZS", # Merchandise trade (% of GDP)
  "EG.CFT.ACCS.ZS", # Access to clean fuels & technologies for cooking (%)
  "FP.CPI.TOTL.ZG", # Inflation, consumer prices (annual %)
  "SM.POP.NETM", # Net migration
  "NE.TRD.GNFS.ZS" # Trade (% of GDP)
)

long_dat <- long_dat %>% dplyr::filter(Series.Code %in% target_codes)

# -------------------------------
# 1) Wide format by series code
# -------------------------------
wide_dat <- long_dat %>%
  tidyr::pivot_wider(
    id_cols = c(Country.Name, Country.Code, year), # country-year IDs
    names_from = Series.Code, # columns by series code
    values_from = value,
    values_fill = NA_real_,
    values_fn = list(
      value = function(.x) {
        v <- suppressWarnings(as.numeric(.x))
        v <- v[!is.na(v)]
        if (length(v) == 0) NA_real_ else v[1] # first non-NA if duplicated
      }
    )
  ) %>%
  dplyr::rename(country = Country.Name, country_code = Country.Code) %>%
  dplyr::relocate(country, country_code, year)

# Clean/short column names used internally
rename_map <- c(
  "NY.GDP.PCAP.KD" = "gdp_pc_const", # per-capita GDP (constant 2015 US$)
  "NY.GDP.MKTP.KD" = "gdp_const", # total GDP (constant 2015 US$)
  "SP.POP.TOTL" = "population",
  "NE.EXP.GNFS.ZS" = "exp_gdp", # Exports (% of GDP)
  "TG.VAL.TOTL.GD.ZS" = "merch_trade_gdp", # Merchandise trade (% of GDP)
  "EG.CFT.ACCS.ZS" = "clean_cooking_access",
  "FP.CPI.TOTL.ZG" = "inflation_cpi", # annual %
  "SM.POP.NETM" = "net_migration",
  "NE.TRD.GNFS.ZS" = "trade_gdp" # Trade (% of GDP)
)

wide_dat <- wide_dat %>%
  dplyr::rename_with(
    ~ rename_map[.],
    .cols = dplyr::any_of(names(rename_map))
  ) %>%
  dplyr::mutate(country_id = dplyr::dense_rank(country)) %>%
  dplyr::relocate(country_id, .before = country)

# Ensure all expected columns exist and are numeric
expected_vars <- c(
  "gdp_pc_const",
  "gdp_const",
  "population",
  "exp_gdp",
  "merch_trade_gdp",
  "trade_gdp",
  "clean_cooking_access",
  "inflation_cpi",
  "net_migration"
)

missing_vars <- setdiff(expected_vars, names(wide_dat))
for (v in missing_vars) {
  wide_dat[[v]] <- NA_real_
}

wide_dat <- wide_dat %>%
  dplyr::mutate(dplyr::across(
    dplyr::all_of(expected_vars),
    ~ suppressWarnings(as.numeric(.x))
  ))

# Compute per-capita GDP when missing (safe fallback, uses MKTP/POP)
wide_dat <- wide_dat %>%
  dplyr::mutate(
    gdp_pc_const = dplyr::coalesce(
      gdp_pc_const,
      dplyr::if_else(
        !is.na(gdp_const) & !is.na(population) & population > 0,
        gdp_const / population,
        NA_real_
      )
    )
  )
country_dict <- wide_dat %>%
  distinct(
    iso3 = country_code,
    country = country
  ) %>%
  filter(!is.na(iso3), !is.na(country))

final_dat <- wide_dat %>%
  dplyr::select(
    country_id,
    year,
    country,
    dplyr::all_of(expected_vars)
  ) %>%
  dplyr::arrange(country_id, year)

# -------------------------------
# 2) (Optional) FRED backfill for inflation
# -------------------------------
fred_key <- Sys.getenv("FRED_API") # standardized name
cc_lu <- raw %>%
  dplyr::distinct(country = Country.Name, country_code = Country.Code)
infl_col <- "inflation_cpi"

# Countries needing inflation backfill
need_fill <- final_dat %>%
  dplyr::left_join(cc_lu, by = "country") %>%
  dplyr::group_by(country_id, country, country_code) %>%
  dplyr::summarise(need = any(is.na(.data[[infl_col]])), .groups = "drop") %>%
  dplyr::filter(need, !is.na(country_code))

fetch_fred_inflation <- function(
  iso3,
  start_year = min(final_dat$year),
  end_year = max(final_dat$year)
) {
  series_id <- paste0("FPCPITOTLZG", iso3)
  url <- "https://api.stlouisfed.org/fred/series/observations"
  out <- tryCatch(
    {
      r <- httr::GET(
        url,
        query = list(
          series_id = series_id,
          api_key = fred_key,
          file_type = "json",
          observation_start = paste0(start_year, "-01-01"),
          observation_end = paste0(end_year, "-12-31"),
          frequency = "a"
        )
      )
      httr::stop_for_status(r)
      jsonlite::fromJSON(httr::content(r, as = "text", encoding = "UTF-8"))
    },
    error = function(e) NULL
  )
  if (is.null(out) || is.null(out$observations)) {
    return(tibble(
      country_code = iso3,
      year = integer(0),
      fred_infl = numeric(0)
    ))
  }
  tibble(
    country_code = iso3,
    year = as.integer(substr(out$observations$date, 1, 4)),
    fred_infl = suppressWarnings(as.numeric(out$observations$value))
  ) %>%
    dplyr::filter(!is.na(fred_infl))
}

if (nrow(need_fill) > 0 && nzchar(fred_key)) {
  fred_df <- purrr::map_dfr(need_fill$country_code, fetch_fred_inflation)
  if (nrow(fred_df) > 0) {
    final_dat <- final_dat %>%
      dplyr::left_join(cc_lu, by = "country") %>%
      dplyr::left_join(fred_df, by = c("country_code", "year")) %>%
      dplyr::mutate(
        inflation_cpi = dplyr::coalesce(inflation_cpi, fred_infl)
      ) %>%
      dplyr::select(-country_code, -fred_infl)
  }
} # else: skip backfill (no missing or no API key)

# -------------------------------
# 3) Scale percentages to proportions
# -------------------------------
final_dat <- final_dat %>%
  dplyr::mutate(
    merch_trade_gdp = merch_trade_gdp / 100,
    exp_gdp = exp_gdp / 100,
    trade_gdp = trade_gdp / 100,
    inflation_cpi = inflation_cpi / 100,
    clean_cooking_access = clean_cooking_access / 100
  )

# -------------------------------
# 4) Sudan / South Sudan rule
# -------------------------------
# 2011–2015: per-capita GDP = (GDP_Sudan + GDP_SouthSudan) / (POP_Sudan + POP_SouthSudan)
# 2000–2010: use Sudan only; South Sudan is not kept anywhere.

sudan_combined_pc <- final_dat %>%
  dplyr::filter(country %in% c("Sudan", "South Sudan"), year >= 2011) %>%
  dplyr::select(country, year, gdp_const, population) %>%
  dplyr::group_by(year) %>%
  dplyr::summarise(
    gdp_sum = if (any(is.na(gdp_const))) NA_real_ else sum(gdp_const),
    pop_sum = if (any(is.na(population))) NA_real_ else sum(population),
    gdp_pc_const_combined = dplyr::if_else(
      !is.na(gdp_sum) & !is.na(pop_sum) & pop_sum > 0,
      gdp_sum / pop_sum,
      NA_real_
    ),
    .groups = "drop"
  )

final_dat <- final_dat %>%
  dplyr::left_join(
    sudan_combined_pc %>% dplyr::select(year, gdp_pc_const_combined),
    by = "year"
  ) %>%
  dplyr::mutate(
    gdp_pc_const = dplyr::if_else(
      country == "Sudan" & year >= 2011,
      gdp_pc_const_combined,
      gdp_pc_const
    )
  ) %>%
  dplyr::select(-gdp_pc_const_combined) %>%
  dplyr::filter(country != "South Sudan") %>% # drop South Sudan entirely
  dplyr::select(-gdp_const, -population) %>% # keep only per-capita after the rule
  dplyr::mutate(country_id = dplyr::dense_rank(country)) %>%
  dplyr::arrange(country_id, year)

# -------------------------------
# 5) Strict drop of any-NA countries
# -------------------------------
target_vars <- c(
  "gdp_pc_const",
  "exp_gdp", # <- added back; remove this line if you intend to exclude it
  "merch_trade_gdp",
  "trade_gdp",
  "clean_cooking_access",
  "inflation_cpi",
  "net_migration"
)

final_dat_strict <- final_dat %>%
  dplyr::group_by(country_id, country) %>%
  dplyr::filter(
    !any(sapply(dplyr::across(dplyr::all_of(target_vars)), anyNA))
  ) %>%
  dplyr::ungroup()

# -------------------------------
# 6) Additionally: WDI-like column names for readability + dictionary
# -------------------------------
# Create a WDI-ish version that respects the original naming style (make.names() style).

clean_to_wdi <- c(
  gdp_pc_const = "GDP.per.capita..constant.2015.US..",
  exp_gdp = "Exports.of.goods.and.services....of.GDP.",
  merch_trade_gdp = "Merchandise.trade....of.GDP.",
  trade_gdp = "Trade....of.GDP.",
  clean_cooking_access = "Access.to.clean.fuels.and.technologies.for.cooking....of.population.",
  inflation_cpi = "Inflation..consumer.prices..annual...",
  net_migration = "Net.migration"
  # (gdp_const, population were dropped by design)
)

present2 <- intersect(names(clean_to_wdi), names(final_dat))

final_dat <- final_dat %>%
  dplyr::rename_with(
    ~ clean_to_wdi[.x],
    .cols = dplyr::all_of(present2)
  )

final_dat_strict <- final_dat_strict %>%
  dplyr::rename_with(
    ~ clean_to_wdi[.x],
    .cols = dplyr::all_of(present2)
  )

# Column dictionary for the replication package
col_dict <- tibble::tibble(
  clean_name = present2,
  wdi_like = unname(clean_to_wdi[present2]),
  note = c(
    gdp_pc_const = "Per-capita GDP (constant 2015 USD)",
    exp_gdp = "Exports of goods & services (% of GDP, proportion after scaling)",
    merch_trade_gdp = "Merchandise trade (% of GDP, proportion)",
    trade_gdp = "Trade (% of GDP, proportion)",
    clean_cooking_access = "Access to clean fuels/tech for cooking (% of population, proportion)",
    inflation_cpi = "Inflation, CPI (annual, proportion)",
    net_migration = "Net migration (persons)"
  )[present2]
)

readr::write_csv(final_dat, "data/processed/Sudan/africa_panel.csv")
readr::write_csv(final_dat_strict, "data/processed/Sudan/africa_panel_drop.csv")
readr::write_csv(col_dict, "data/processed/Sudan/column_dictionary.csv")

# ============================================================
# DOT -> trade matrix (no row-normalization, no single-side fill)
#   1) build annual bilateral amounts A_ij(t) using ONLY years
#      where BOTH directions (i->j and j->i) are observed
#   2) restrict countries to those in final_dat$country
#   3) extract Sudan row vector from the averaged (unnormalized) matrix
#   4) output one matrix (A_avg) and one vector (Sudan row)
# Source: International Monetary Fund
# ============================================================

library(tidyverse)

# --- Parameters
years_pre <- 2000:2010
sudan_name <- "Sudan"

# --- Load DOT raw (replace path)
dot_raw <- read.csv("data/raw/IMF_trade.csv", stringsAsFactors = FALSE)

dot_iso <- dot_raw %>%
  mutate(
    reporter_iso3 = countrycode::countrycode(
      COUNTRY,
      origin = "country.name",
      destination = "iso3c",
      warn = FALSE
    ),
    partner_iso3 = countrycode::countrycode(
      COUNTERPART_COUNTRY,
      origin = "country.name",
      destination = "iso3c",
      warn = FALSE
    )
  )

countries_keep_iso3 <- unique(raw$Country.Code)
# --- Scale map (extend if needed)
scale_map <- c("Millions" = 1e6, "Units" = 1)

# --- Clean & parse
dot_clean <- dot_iso %>%
  filter(
    reporter_iso3 %in% countries_keep_iso3,
    partner_iso3 %in% countries_keep_iso3
  ) %>%
  rename(
    reporter = reporter_iso3,
    partner = partner_iso3,
    year = TIME_PERIOD
  )

dot_clean <- dot_clean %>%
  mutate(
    multiplier = unname(scale_map[trimws(SCALE)]),
    multiplier = ifelse(is.na(multiplier), 1, multiplier),
    value_usd = as.numeric(OBS_VALUE) * multiplier,
    flow = case_when(
      grepl("^Exports", INDICATOR, ignore.case = TRUE) ~ "export",
      grepl("^Imports", INDICATOR, ignore.case = TRUE) ~ "import",
      TRUE ~ NA_character_
    ),
    valuation = case_when(
      grepl("CIF", INDICATOR, ignore.case = TRUE) ~ "CIF",
      grepl("FOB", INDICATOR, ignore.case = TRUE) ~ "FOB",
      TRUE ~ NA_character_
    )
  ) %>%
  mutate(year = as.integer(year)) %>%
  filter(FREQUENCY == "Annual", year %in% years_pre)

# --- Collapse valuation/flow per (i,j,t)
#     imports: prefer CIF > FOB; exports: FOB; then trade_ij = imports_pref + exports_fob
flows_by_val <- dot_clean %>%
  group_by(reporter, partner, year, flow, valuation) %>%
  summarise(value_usd = sum(value_usd, na.rm = TRUE), .groups = "drop")

wide_flows <- flows_by_val %>%
  mutate(col = paste0(flow, "_", ifelse(is.na(valuation), "NA", valuation))) %>%
  select(reporter, partner, year, col, value_usd) %>%
  pivot_wider(names_from = col, values_from = value_usd, values_fill = 0)

trade_ij_year <- wide_flows %>%
  mutate(
    imports_pref = ifelse(
      !is.na(import_CIF) & import_CIF > 0,
      import_CIF,
      ifelse(!is.na(import_FOB) & import_FOB > 0, import_FOB, NA_real_)
    ),
    exports_fob = ifelse(
      !is.na(export_FOB) & export_FOB > 0,
      export_FOB,
      NA_real_
    ),
    trade_ij = ifelse(
      is.na(imports_pref) & is.na(exports_fob),
      NA_real_,
      coalesce(imports_pref, 0) + coalesce(exports_fob, 0)
    )
  ) %>%
  select(reporter, partner, year, trade_ij)

# --- Require BOTH directions in a year; drop single-sided observations
dir_ij <- trade_ij_year %>% rename(i = reporter, j = partner, v_ij = trade_ij)
dir_ji <- trade_ij_year %>% rename(j = reporter, i = partner, v_ji = trade_ij)

sym_both <- inner_join(dir_ij, dir_ji, by = c("i", "j", "year")) %>%
  mutate(A_ij = (v_ij + v_ji) / 2) %>%
  select(i, j, year, A_ij) %>%
  filter(!is.na(A_ij), i != j)

# --- Build annual matrices with NA (no filling)
make_matrix <- function(df_y, countries) {
  M <- matrix(
    NA_real_,
    nrow = length(countries),
    ncol = length(countries),
    dimnames = list(countries, countries)
  )
  if (nrow(df_y) > 0) {
    for (k in seq_len(nrow(df_y))) {
      i <- df_y$i[k]
      j <- df_y$j[k]
      a <- df_y$A_ij[k]
      if (!is.na(a) && i %in% countries && j %in% countries && i != j) {
        M[i, j] <- a
        M[j, i] <- a
      }
    }
  }
  diag(M) <- 0
  M
}

filtered_country <- final_dat_strict %>% select(country) %>% distinct()
countries_keep_iso3 <- country_dict %>%
  filter(country %in% filtered_country$country)
countries_keep_iso3 <- countries_keep_iso3$iso3

A_year_list <- sym_both %>%
  dplyr::group_by(year) %>%
  dplyr::group_map(~ make_matrix(.x, countries_keep_iso3)) %>%
  setNames(sort(unique(sym_both$year)))

if (length(A_year_list) == 0) {
  stop("No bilateral pairs with both directions observed in the pre period.")
}

# --- Pre-intervention average WITHOUT filling:
#     average over years where A_ij(t) is observed; if never observed -> NA
A_sum <- Reduce(
  "+",
  lapply(A_year_list, function(M) {
    M[is.na(M)] <- 0
    M
  })
)
A_count <- Reduce(
  "+",
  lapply(A_year_list, function(M) {
    (!is.na(M)) * 1
  })
)
A_avg <- A_sum / ifelse(A_count > 0, A_count, NA_real_)
dimnames(A_avg) <- list(countries_keep_iso3, countries_keep_iso3)
diag(A_avg) <- 0

# --- Zero-fill ONLY pairs never observed in any year (keep observed means)
A_avg0 <- A_avg
A_avg0[is.na(A_avg0)] <- 0
diag(A_avg0) <- 0

# --- Sudan row vector (unnormalized, after zero-fill)
sudan_name <- "SDN"
w_sudan <- if (sudan_name %in% rownames(A_avg0)) {
  A_avg0[sudan_name, , drop = TRUE]
} else {
  NULL
}

# --- Drop Sudan from the matrix (both row and column) and from the vector
countries <- rownames(A_avg0)
sudan_idx <- which(countries == sudan_name)

if (length(sudan_idx) > 0) {
  # 1) remove Sudan row/column from the matrix
  A_avg0 <- A_avg0[-sudan_idx, -sudan_idx, drop = FALSE]

  # 2) remove Sudan element from the Sudan vector (to align with the matrix columns)
  sudan_vec_final <- if (!is.null(w_sudan)) w_sudan[-sudan_idx] else NULL

  # names stay aligned by subsetting; set explicitly for safety
  rownames(A_avg0) <- colnames(A_avg0) <- countries[-sudan_idx]
  if (!is.null(sudan_vec_final)) names(sudan_vec_final) <- countries[-sudan_idx]
} else {
  # Sudan already excluded from the matrix; keep vector as is (may be NULL)
  sudan_vec_final <- w_sudan
}

iso_to_name <- setNames(country_dict$country, country_dict$iso3)
country_iso <- colnames(A_avg0)
country_names <- unname(iso_to_name[country_iso])
dimnames(A_avg0) <- list(country_names, country_names)
# --- Write outputs
write.csv(A_avg0, "data/processed/Sudan/weight_mat.csv", row.names = TRUE)
if (!is.null(sudan_vec_final)) {
  write.csv(
    data.frame(
      country = colnames(A_avg0),
      amount = as.numeric(sudan_vec_final)
    ),
    "data/processed/Sudan/weight_vec.csv",
    row.names = FALSE
  )
}

rm(list = ls())
panel_df <- read_csv("data/processed/Sudan/africa_panel_drop.csv")
w <- read_csv("data/processed/Sudan/weight_vec.csv") %>% arrange(country)
W <- read_csv("data/processed/Sudan/weight_mat.csv") %>% arrange("...1")

countries_panel <- panel_df %>%
  distinct(country_id, country) %>%
  arrange(country_id) %>%
  pull(country)

countries_W <- W$...1
countries_w <- w$country

common_countries <- Reduce(
  intersect,
  list(countries_panel, countries_W, countries_w)
)

panel_df <- panel_df %>%
  filter(country %in% common_countries | country == "Sudan")

order_countries <- countries_panel[countries_panel %in% common_countries]

W_mat <- W %>%
  column_to_rownames("...1") %>%
  as.matrix()

W_mat <- W_mat[order_countries, order_countries, drop = FALSE]

W <- W_mat %>%
  as.data.frame() %>%
  rownames_to_column("country") %>%
  as_tibble()

w <- w %>%
  filter(country %in% order_countries) %>%
  mutate(country = factor(country, levels = order_countries)) %>%
  arrange(country)

sudan_secession <- list(
  panel_df = panel_df,
  w = w,
  W = W
)

usethis::use_data(sudan_secession, overwrite = TRUE)
