# Variable Dictionary

This file documents all variables included in the bundled datasets used by the replication scripts.

## `data/california_smoking.rda`

Object contains three data frames: `panel`, `w`, `W`.

### `panel`
- `state`: US state name.
- `year`: Calendar year.
- `cigsale`: Per-capita cigarette sales outcome.
- `retprice`: Retail cigarette price covariate.
- `state_id`: Numeric state identifier used in preprocessing.

### `w`
- `state`: Control-state name.
- `w`: Spillover exposure weight vector for the treated unit.

### `W`
- First column is the row unit/state identifier.
- Remaining columns define the row-normalized spatial adjacency matrix among controls.

## `data/sudan_secession.rda`

Object contains three data frames: `panel_df`, `w`, `W`.

### `panel_df`
- `country_id`: Numeric country identifier.
- `year`: Calendar year.
- `country`: Country name.
- `GDP.per.capita..constant.2015.US..`: GDP per capita (constant 2015 USD), outcome variable.
- `Exports.of.goods.and.services....of.GDP.`: Exports (% of GDP).
- `Merchandise.trade....of.GDP.`: Merchandise trade (% of GDP).
- `Trade....of.GDP.`: Total trade (% of GDP).
- `Access.to.clean.fuels.and.technologies.for.cooking....of.population.`: Access to clean cooking fuels (% population).
- `Inflation..consumer.prices..annual...`: CPI inflation (annual %).
- `Net.migration`: Net migration.

### `w`
- `country`: Control-country name.
- `w`: Spillover exposure weight vector for treated Sudan.

### `W`
- First column is the row unit/country identifier.
- Remaining columns define the row-normalized spatial/trade-weight matrix among controls.

## Non-proprietary copies

CSV copies of each dataset are provided under `data/nonproprietary/`:
- `california_panel.csv`
- `california_w_vector.csv`
- `california_W_matrix.csv`
- `sudan_panel.csv`
- `sudan_w_vector.csv`
- `sudan_W_matrix.csv`
