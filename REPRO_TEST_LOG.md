# Reproducibility Test Log

## Environment
- Date: 2026-02-10
- Command mode: smoke
- R: 4.2.3
- Platform: Darwin 25.3.0 arm64

## Command executed

```bash
SCSPILL_MODE=smoke SCSPILL_CLEAN_OUTPUT=true Rscript code/99_run_all.R smoke all
```

## Result
- Status: success (exit code 0)
- Logs:
  - `output/logs/20260210-134922_01_california_main_smoke.log`
  - `output/logs/20260210-134922_02_sudan_main_smoke.log`
  - `output/logs/20260210-134922_03_simulation_main_smoke.log`
  - `output/logs/20260210-134922_04_geweke_main_smoke.log`

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
- No fatal errors and no execution halts.

## Submission-boundary note
- Figures/tables generated during this smoke run are expected reproducible artifacts and are **not bundled** for submission.
- After test completion, `output/figures/` and `output/tables/` were reset to placeholder-only state (`.gitkeep`) to keep the submission boundary explicit.
