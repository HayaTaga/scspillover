# ECTJ Reproducibility Requirements Audit

Source documents treated as authoritative specification:
- `ECTJ_rep_1.pdf` (Before submission / policy)
- `ECTJ_rep2.pdf` (Best practices)
- `ECTJ_rep3.pdf` (Prepare package / README details)

## Step 1. Extracted requirements

### Mandatory

| ID | Requirement (short) | PDF wording/section | Verifiable acceptance criterion |
|---|---|---|---|
| M1 | Complete replication package is required for publication | Policy text: “We will only publish papers if they are accompanied by a complete replication package.” | Package contains data + code + README and can reproduce paper outputs. |
| M2 | README must exist and describe package + software versions | Replication Policy + Package §2/§3 | `README_submit.md` and `README.pdf` document structure, software versions, packages. |
| M3 | Include data used in analysis with transparent variable documentation | Policy + Package §2 | Datasets included under `data/`; variable dictionary present (`data/VARIABLE_DICTIONARY.md`). |
| M4 | Include all code to reproduce tables/figures/supplement/appendix | Policy + Package §2 | Scripts in `code/` generate all mapped outputs listed in README mapping table. |
| M5 | For Monte Carlo/simulation, set seeds for exact reproducibility | Policy + Package §2 | Scripts contain explicit `set.seed(...)` / fixed seeds. |
| M6 | Include non-proprietary copy when data are in proprietary format | Policy + Package §2 | CSV copies included under `data/nonproprietary/`. |
| M7 | README must include data availability statement | Package §3 item 2 | Section present in `README_submit.md`. |
| M8 | README must include rights statements (both certifications) | Package §3 item 3 | Exact certification statements included in `README_submit.md`. |
| M9 | README must include run instructions + where outputs are stored | Package §3 items 4-5 | Master command + output directories documented in `README_submit.md`. |
| M10 | README must include software/OS requirements + package install information | Package §3 items 6-7 | Requirements table and dependency instructions in `README_submit.md`; `code/00_setup.R`. |
| M11 | README must include expected run time (+ hardware where relevant) | Package §3 item 8 | Runtime estimates by script are documented. |
| M12 | Data citations in paper references + copy in README | Package §3 item 9 and §4 | Data citation section included in `README_submit.md`. |

### Recommended (best practices)

| ID | Requirement (short) | PDF wording/section | Verifiable acceptance criterion |
|---|---|---|---|
| R1 | Separate code/data/output directories | Best Practices (Folder Structure, Code, Output) | Distinct `code/`, `data/`, `output/`. |
| R2 | Use scripts (not manual commands), with master file | Best Practices (Code) | `code/99_run_all.R` calls subsidiary scripts. |
| R3 | Run all-at-once design, paths set once/relative | Best Practices (Code) | Single entry execution from repository root. |
| R4 | Log outputs to disk | Best Practices (Output) | Execution logs generated under `output/logs/`. |
| R5 | Include table mapping paper outputs to scripts/files | Best Practices (Output) | Mapping table in `README_submit.md`. |
| R6 | Test from clean folder with expected output deleted | Best Practices (Output) | `code/99_run_all.R` optionally cleans `output/` before run. |
| R7 | Keep raw data intact and track versions | Best Practices (Data) | Source/provenance documented; exact bundled extracts versioned in repo. |

### Optional

| ID | Requirement (short) | PDF wording/section | Verifiable acceptance criterion |
|---|---|---|---|
| O1 | Provide makefile | Policy: “Ideally, but not necessarily, include a makefile” | `Makefile` included. |
| O2 | Paper source depends directly on generated outputs | Best Practices (Output) | Documented aspiration; not mandatory for check pass. |

## Step 2. Repository audit before refactor

| ID | Status before refactor | Evidence |
|---|---|---|
| M1 | Partially satisfied | Code/data present, but reproducibility depended on large `vignettes/*.Rmd` and inconsistent output mapping. |
| M2 | Not satisfied | No `README.pdf`; only `README.md` and an outdated `README_submit.md`. |
| M3 | Partially satisfied | `.rda` datasets existed, but no standalone variable dictionary. |
| M4 | Partially satisfied | Code existed, but spread in Rmd notebooks; no robust master execution script. |
| M5 | Satisfied | Fixed seeds present in scripts (vignettes and R functions). |
| M6 | Not satisfied | No CSV non-proprietary copies for bundled `.rda` objects. |
| M7 | Partially satisfied | Data availability existed but not fully aligned with ECTJ README requirements. |
| M8 | Not satisfied | Required rights certifications were missing verbatim. |
| M9 | Partially satisfied | Instructions existed, but referenced non-existing paths (`Rmarkdown/`) and inconsistent outputs. |
| M10 | Partially satisfied | Dependencies partially documented; no preflight checker script. |
| M11 | Partially satisfied | Runtime estimates existed but without clear mode separation and execution entrypoint. |
| M12 | Partially satisfied | Some references existed; dedicated ECTJ-style data citation section incomplete. |
| R1 | Partially satisfied | Code/data existed, but outputs were mixed under `inst/` and pre-generated artifacts tracked. |
| R2 | Not satisfied | No master script orchestrating full workflow. |
| R3 | Partially satisfied | Relative paths used in places, but run flow relied on manual notebook rendering. |
| R4 | Not satisfied | No structured logs directory/workflow. |
| R5 | Partially satisfied | Mapping table existed but had wrong script paths and mismatched output filenames. |
| R6 | Not satisfied | No automated clean-run mechanism. |
| R7 | Partially satisfied | Provenance narrative existed but exact dictionary/format export missing. |
| O1 | Not satisfied | No Makefile in repository root. |

## Step 3. Minimal sufficient refactoring plan (implemented)

1. Move substantive computation into plain `.R` scripts under `code/` and keep `vignettes/` as thin wrappers.
2. Add master execution script and optional Makefile.
3. Separate generated artifacts into `output/` and stop relying on `inst/` outputs.
4. Add package preflight setup script and run modes (`full`/`smoke`) for reproducibility checks.
5. Add non-proprietary CSV exports and variable dictionary.
6. Rewrite submission README with strict output mapping and explicit assumptions.
7. Remove compiled binaries (`src/*.o`, `src/*.so`) from reliance path.

## Step 4. Final consistency check (submission-hardening)

### Mandatory requirements (M1-M12)

| ID | Final status | Evidence |
|---|---|---|
| M1 | satisfied | `code/99_run_all.R` + `RESULTS_MAPPING.csv` define complete reproducible workflow. |
| M2 | satisfied | `README_submit.md` and `README.pdf` exist; software/dependency versions documented with `DEPENDENCY_LOCK.csv`. |
| M3 | satisfied | Input data in `data/` including `data/raw/`, bundled `.rda`, and variable documentation in `data/VARIABLE_DICTIONARY.md`. |
| M4 | satisfied | Reproduction scripts in `code/01_*.R` to `code/04_*.R`. |
| M5 | satisfied | Fixed seed initialization in scripts and seed contracts recorded in logs. |
| M6 | satisfied | Non-proprietary dataset copies in `data/nonproprietary/`. |
| M7 | satisfied | Data availability statement in `README_submit.md` Section 3. |
| M8 | satisfied | Required rights statements in `README_submit.md` Section 4. |
| M9 | satisfied | Run instructions, mode contract, output mapping, and submission boundary in `README_submit.md` Sections 6-9. |
| M10 | satisfied | Software/package requirements and lock handling in `README_submit.md` Section 5 and `code/00_setup.R`. |
| M11 | satisfied | Runtime/hardware guidance in `README_submit.md` Section 6 (`full` mode). |
| M12 | satisfied | Data citations in `README_submit.md` Section 11. |

### Recommended requirements (R1-R6 as implemented)

| ID | Final status | Evidence |
|---|---|---|
| R1 | satisfied | `code/`, `data/`, `output/` are separated. |
| R2 | satisfied | Master script `code/99_run_all.R` orchestrates sub-scripts. |
| R3 | satisfied | Single-entry scripted run from repository root, relative paths only. |
| R4 | satisfied | Disk logs written to `output/logs/` with environment metadata. |
| R5 | satisfied | Output-to-script mapping in `RESULTS_MAPPING.csv` and documented in `README_submit.md`. |
| R6 | satisfied | Clean-output behavior in `code/99_run_all.R` via `SCSPILL_CLEAN_OUTPUT`. |
