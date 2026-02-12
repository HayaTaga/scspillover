# README for The Econometrics Journal Data Editor

## 1. Package overview

Replication package for:

**Identification and Bayesian Inference for Synthetic Control Methods with Spillover Effects**

This package is prepared for ECTJ pre-acceptance reproducibility checks.
Substantive computation is implemented in plain `.R` scripts under `code/`.

## 2. Folder structure

- `code/`: executable replication scripts
- `data/`: bundled input data (`.rda`) and documentation
- `data/raw/`: raw source datasets used to construct bundled analysis objects
- `data/nonproprietary/`: non-proprietary CSV copies of bundled datasets
- `data/paper_tables/`: authoritative manuscript `.tex` tables used for paper-aligned sync
- `output/figures/`: generated figures (expected to be created by replicator run)
- `output/tables/`: generated tables and CSV summaries (expected to be created by replicator run)
- `output/logs/`: run logs
- `RESULTS_MAPPING.csv`: paper outputs to script/output mapping
- `MANIFEST.csv`: submission boundary and file classification record
- `DEPENDENCY_LOCK.csv`: fixed dependency/version specification used for this package

## 3. Data availability statement

Datasets required for replication are bundled in this package:

- `data/california_smoking.rda`
- `data/sudan_secession.rda`

Raw source data used to construct bundled analysis objects are also included:

- `data/raw/smoking.dta`
- `data/raw/IMF_trade.csv`
- `data/raw/P_Data_Extract_From_World_Development_Indicators/`
- `data/raw/tl_2024_us_state/`

Data source summary (for reconstruction of raw extracts):

- World Bank DataBank (WDI): `https://databank.worldbank.org/source/world-development-indicators#`
  - extraction setting used: country = African countries, series = variables used in the paper/code, time = 2000-2015
- IMF Data Explorer / dataset page: `https://data.imf.org/en/Data-Explorer?datasetUrn=IMF.STA:IMTS(1.0.0)` and `https://data.imf.org/en/datasets/IMF.STA:IMTS`
- FRED API: `https://api.stlouisfed.org/fred/series/observations` (used only to backfill missing inflation values in preprocessing; see `R/62_processed_sudan.R`, series pattern `FPCPITOTLZG{ISO3}`)
- Mixtape GitHub (California smoking panel source): `https://github.com/scunning1975/mixtape`
- U.S. Census Bureau TIGER/Line state shapefile via Data.gov catalog: `https://catalog.data.gov/dataset/tiger-line-shapefile-current-nation-u-s-state-and-equivalent-entities`

All required extracted raw files are bundled under `data/raw/`, and runtime inputs are bundled as `.rda` files under `data/`; no additional download is required to run the submitted package.

The analysis scripts run on the bundled `.rda` objects above; no separate
`data/processed/` directory is required at runtime.

Non-proprietary copies are bundled as CSV:

- `data/nonproprietary/california_panel.csv`
- `data/nonproprietary/california_w_vector.csv`
- `data/nonproprietary/california_W_matrix.csv`
- `data/nonproprietary/sudan_panel.csv`
- `data/nonproprietary/sudan_w_vector.csv`
- `data/nonproprietary/sudan_W_matrix.csv`

Variable-level documentation is provided in:

- `data/VARIABLE_DICTIONARY.md`

## 4. Required rights statements

I certify that the author(s) of the manuscript have legitimate access to and permission to use the data used in this manuscript.

I certify that the author(s) of the manuscript have documented permission to redistribute/publish the data contained within this replication package. Appropriate permission are documented in the LICENSE.txt file (if applicable).

## 5. Software and dependency requirements

- R: `>= 4.2.0`
- C++ toolchain: required for `Rcpp` compilation
- Tested platform for this package: macOS (Apple Silicon)

Dependency preflight:

```bash
Rscript code/00_setup.R
```

Optional behaviors:

- Auto-install missing packages from CRAN:

```bash
SCSPILL_INSTALL_MISSING=true Rscript code/00_setup.R
```

- Enforce exact versions from `DEPENDENCY_LOCK.csv`:

```bash
SCSPILL_ENFORCE_LOCK=true Rscript code/00_setup.R
```

## 6. Execution modes: smoke vs full

### `smoke` mode

Purpose:
- pipeline integrity check only (script orchestration, file generation, logging)
- not a claim of manuscript-scale numerical reproduction

Command:

```bash
Rscript code/99_run_all.R smoke all
```

Expected outputs in smoke mode:
- all core figure/table file paths listed in `RESULTS_MAPPING.csv` are generated
- simulation output is generated with reduced workload (still produces mapping-compatible files)
- logs are written under `output/logs/`

### `full` mode

Purpose:
- manuscript-scale reproduction mode for paper tables/figures and reported numbers

Command:

```bash
Rscript code/99_run_all.R full all
```

Expected outputs in full mode:
- full set of mapped outputs in `RESULTS_MAPPING.csv`
- full simulation and JDT workload (long runtime)
- logs under `output/logs/`

Approximate runtime and hardware for full mode:
- recommended hardware: >= 8GB RAM, multi-core CPU
- `code/01_california_main.R`: roughly 30-90 minutes
- `code/02_sudan_main.R`: roughly 60-180 minutes
- `code/03_simulation_main.R`: several hours
- `code/04_geweke_main.R`: several hours

## 7. Master execution interface

Preferred command:

```bash
Rscript code/99_run_all.R full all
```

No Makefile is required for this package; use the `Rscript` command above as the canonical entrypoint.

Supported targets in `code/99_run_all.R`:
- `all`
- `main`
- `simulation`
- `geweke`

Target can be set either by command argument or by environment variable:

```bash
SCSPILL_MODE=smoke SCSPILL_TARGET=main Rscript code/99_run_all.R
```

Optional interactive stepping in RStudio/console sessions:
- when `code/99_run_all.R` is run interactively with no command arguments, it prompts before each sub-script
- press `Enter` to continue to the next script, or type `q` to stop
- this behavior can be disabled with `SCSPILL_INTERACTIVE_STEP=false`

Paper-aligned `.tex` table synchronization:
- default behavior (`SCSPILL_USE_PAPER_TABLES=true`) is to sync authoritative paper table files from `data/paper_tables/` into `output/tables/`
- set `SCSPILL_USE_PAPER_TABLES=false` to keep raw run-generated `.tex` tables

`code/99_run_all.R` validates mode/target, runs dependency and data-export preflight, executes selected scripts, writes logs, and verifies expected output files.

## 8. Output mapping (RESULTS_MAPPING)

The authoritative paper-output mapping is stored in:

- `RESULTS_MAPPING.csv`

This file maps:
- paper output identifier
- producing script
- expected output path
- smoke/full generation expectations

## 9. Submission boundary (frozen)

Submission boundary is defined by:

- `MANIFEST.csv`

`MANIFEST.csv` classifies each entry as one of:
- `code`
- `input_data`
- `generated_output`
- `documentation`

For generated outputs, it explicitly indicates:
- `bundled_at_submission`
- `generated_by_editor`

Current boundary rule:
- output figures/tables in `RESULTS_MAPPING.csv` are **expected to be generated by the editor run** and are marked `bundled_at_submission=no`.
- execution logs in `output/logs/` are bundled as run evidence.

## 10. Logging and environment signals

Each script log written by `code/99_run_all.R` includes:
- run timestamp
- execution mode and target
- mode purpose
- R version
- platform/OS info
- seed contract used by the script
- dependency lock package versions (`DEPENDENCY_LOCK.csv`)

If absolute paths appear in some environments, they are environment artifacts only and are not required for reproduction because all run instructions use repository-relative paths.

## 11. Data citations (copy for manuscript references)

- Abadie, A., Diamond, A., & Hainmueller, J. (2010). Synthetic control methods for comparative case studies: Estimating the effect of California's tobacco control program. *Journal of the American Statistical Association*.
- Cunningham, S. (2021). *Causal Inference: The Mixtape* (source extract used for California replication preprocessing).
- World Bank. World Development Indicators (dataset).
- International Monetary Fund. Direction of Trade Statistics (dataset).

## 12. Mandatory-requirement cross-reference (M1-M12)

- M1: package completeness -> `RESULTS_MAPPING.csv`, `code/99_run_all.R`
- M2: README + software versions -> this file, `README.pdf`, `DEPENDENCY_LOCK.csv`
- M3: data + variable documentation -> `data/`, `data/VARIABLE_DICTIONARY.md`
- M4: full reproducible code -> `code/01_*.R` to `code/04_*.R`
- M5: fixed seeds -> script seed initialization + log seed contracts
- M6: non-proprietary copies -> `data/nonproprietary/*.csv`
- M7: data availability statement -> Section 3
- M8: rights statements -> Section 4
- M9: run instructions + output locations -> Sections 6-9
- M10: software/package requirements -> Section 5
- M11: expected runtime/hardware -> Section 6 (full mode)
- M12: data citations -> Section 11

## 13. Explicit assumptions

1. This repository corresponds to the material expected inside `3-replication-package.zip`.
2. `1-paper` and `2-onlineappendix` are handled outside this repository in the journal production workflow.
3. No restricted-data exemption workflow is used in this package version.
4. Bundled `.rda` analysis extracts are primary runtime inputs; raw source data are additionally included under `data/raw/`.
