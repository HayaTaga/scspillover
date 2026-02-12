# Reproducibility Test Log

## Environment
- Date: 2026-02-11
- Command mode: smoke
- R: 4.5.0
- Platform: Darwin 25.3.0 arm64

## Command executed

```bash
SCSPILL_MODE=smoke SCSPILL_TARGET=all SCSPILL_CLEAN_OUTPUT=true Rscript code/99_run_all.R smoke all
```

## Result
- Status: success (exit code 0)
- Logs:
  - `output/logs/20260211-160357_01_california_main_smoke.log`
  - `output/logs/20260211-160357_02_sudan_main_smoke.log`
  - `output/logs/20260211-160357_03_simulation_main_smoke.log`
  - `output/logs/20260211-160357_04_geweke_main_smoke.log`

## Log coverage
Each log includes:
- run timestamp
- mode and mode purpose
- R version
- platform/OS
- dependency lock versions
- seed contract and script-level seed initialization lines

## Warnings
- Non-fatal warnings observed in main scripts:
  - `row names were found from a short variable and have been discarded`
  - `Removed 1 row containing missing values or values outside the scale range (geom_ribbon())`
- No fatal errors and no execution halts.

## Submission-boundary note
- Figures/tables generated during this smoke run are expected reproducible artifacts and are **not bundled** for submission.
- If preparing the final archive, remove expected generated files in `output/figures/` and `output/tables/` and keep only folders/placeholders plus logs.
