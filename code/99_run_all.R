args <- commandArgs(trailingOnly = TRUE)
mode <- if (length(args) >= 1) tolower(args[[1]]) else tolower(Sys.getenv("SCSPILL_MODE", "full"))
target <- if (length(args) >= 2) tolower(args[[2]]) else tolower(Sys.getenv("SCSPILL_TARGET", "all"))

if (!mode %in% c("full", "smoke")) {
  stop("Mode must be 'full' or 'smoke'.")
}
if (!target %in% c("all", "main", "simulation", "geweke")) {
  stop("Target must be one of: all, main, simulation, geweke.")
}

interactive_step <- interactive() &&
  length(args) == 0 &&
  identical(tolower(Sys.getenv("SCSPILL_INTERACTIVE_STEP", "true")), "true")

mode_purpose <- switch(
  mode,
  smoke = "Pipeline integrity check only (not a claim of final manuscript-scale numerical reproduction).",
  full = "Manuscript-scale reproduction mode for paper tables/figures and reported numbers."
)
message(sprintf("Execution mode: %s", mode))
message(sprintf("Mode purpose: %s", mode_purpose))

setup_out <- system2(
  "Rscript",
  c("code/00_setup.R"),
  env = c(
    paste0("SCSPILL_MODE=", mode),
    paste0("SCSPILL_INSTALL_MISSING=", Sys.getenv("SCSPILL_INSTALL_MISSING", "false")),
    paste0("SCSPILL_ENFORCE_LOCK=", Sys.getenv("SCSPILL_ENFORCE_LOCK", "false"))
  ),
  stdout = TRUE,
  stderr = TRUE
)
if (!is.null(attr(setup_out, "status")) && attr(setup_out, "status") != 0) {
  cat(setup_out, sep = "\n")
  stop("Dependency setup failed.")
}
cat(setup_out, sep = "\n")

export_out <- system2(
  "Rscript",
  c("code/05_export_nonproprietary_data.R"),
  stdout = TRUE,
  stderr = TRUE
)
if (!is.null(attr(export_out, "status")) && attr(export_out, "status") != 0) {
  cat(export_out, sep = "\n")
  stop("Non-proprietary data export failed.")
}
cat(export_out, sep = "\n")

scripts <- switch(
  target,
  all = c("code/01_california_main.R", "code/02_sudan_main.R", "code/03_simulation_main.R", "code/04_geweke_main.R"),
  main = c("code/01_california_main.R", "code/02_sudan_main.R"),
  simulation = c("code/03_simulation_main.R"),
  geweke = c("code/04_geweke_main.R")
)

script_seed_contract <- c(
  "code/01_california_main.R" = "global set.seed(20251022) and sc_spillover seed=20251022",
  "code/02_sudan_main.R" = "global set.seed(20251022) and sc_spillover seed=20251022",
  "code/03_simulation_main.R" = "global set.seed(20251030) and deterministic per-scenario seed vectors",
  "code/04_geweke_main.R" = "global set.seed(20251030)"
)

read_lock <- function(path = "DEPENDENCY_LOCK.csv") {
  if (!file.exists(path)) {
    return(data.frame(type = character(0), name = character(0), version = character(0), stringsAsFactors = FALSE))
  }
  utils::read.csv(path, stringsAsFactors = FALSE)
}

lock_df <- read_lock("DEPENDENCY_LOCK.csv")

build_log_header <- function(script, mode, target, lock_df, seed_contract) {
  runtime_stamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S %z")
  sys <- Sys.info()
  pkg_rows <- subset(lock_df, type == "package")
  pkg_lines <- if (nrow(pkg_rows) > 0) {
    paste0("- ", pkg_rows$name, "==", pkg_rows$version)
  } else {
    "- (DEPENDENCY_LOCK.csv not found)"
  }

  c(
    sprintf("run_timestamp: %s", runtime_stamp),
    sprintf("script: %s", script),
    sprintf("mode: %s", mode),
    sprintf("target: %s", target),
    sprintf("mode_purpose: %s", mode_purpose),
    sprintf("r_version: %s", R.version.string),
    sprintf("platform: %s %s %s", sys[["sysname"]], sys[["release"]], sys[["machine"]]),
    sprintf("seed_contract: %s", seed_contract),
    "dependency_lock:",
    pkg_lines,
    "--- script output ---"
  )
}

check_expected_outputs <- function(target) {
  main_outputs <- c(
    "output/figures/fig_ca_diag.png",
    "output/figures/fig_ca_panelA.pdf",
    "output/figures/fig_ca_panelB.pdf",
    "output/figures/fig_ca_panel_spillover.pdf",
    "output/figures/fig_ca_prior_predictive.pdf",
    "output/figures/fig_ca_weight.pdf",
    "output/figures/fig_sudan_diag.png",
    "output/figures/fig_sudan_panelA.pdf",
    "output/figures/fig_sudan_panelB.pdf",
    "output/figures/fig_sudan_panel_spillover.pdf",
    "output/figures/fig_sudan_prior_predictive.pdf",
    "output/figures/fig_sudan_trade.pdf",
    "output/figures/fig_sudan_weight.pdf",
    "output/tables/ca_diag_table.tex",
    "output/tables/ca_sens_table.tex",
    "output/tables/ca_ppa_summary_table.tex",
    "output/tables/sudan_diag_table.tex",
    "output/tables/sudan_sens_table.tex",
    "output/tables/sudan_ppa_summary_table.tex"
  )

  simulation_outputs <- c("output/tables/simulation_results.tex")
  geweke_outputs <- c("output/tables/geweke_jdt_summary.tex", "output/tables/geweke_jdt_summary.csv")

  required <- switch(
    target,
    all = c(main_outputs, simulation_outputs, geweke_outputs),
    main = main_outputs,
    simulation = simulation_outputs,
    geweke = geweke_outputs
  )

  missing <- required[!file.exists(required)]
  if (length(missing) > 0) {
    stop(sprintf("Missing expected output files: %s", paste(missing, collapse = ", ")))
  }

  if (target %in% c("all", "simulation")) {
    mc_files <- list.files("output/tables/mc_result", pattern = "^mc_study_.*[.]csv$", full.names = TRUE)
    if (length(mc_files) < 1) {
      stop("Missing expected Monte Carlo result CSVs in output/tables/mc_result/.")
    }
  }
}

dir.create("output/figures", recursive = TRUE, showWarnings = FALSE)
dir.create("output/tables", recursive = TRUE, showWarnings = FALSE)
dir.create("output/logs", recursive = TRUE, showWarnings = FALSE)

clean_output <- identical(Sys.getenv("SCSPILL_CLEAN_OUTPUT", "true"), "true")
if (clean_output) {
  unlink(list.files("output/figures", full.names = TRUE), force = TRUE)
  unlink(list.files("output/tables", full.names = TRUE), recursive = TRUE, force = TRUE)
  dir.create("output/tables/mc_result", recursive = TRUE, showWarnings = FALSE)
}

timestamp <- format(Sys.time(), "%Y%m%d-%H%M%S")
aborted_by_user <- FALSE
for (script in scripts) {
  if (interactive_step) {
    ans <- tolower(trimws(readline(sprintf(
      "About to run %s [mode=%s, target=%s]. Press Enter to continue or type 'q' to stop: ",
      script,
      mode,
      target
    ))))
    if (ans %in% c("q", "quit")) {
      message(sprintf("Stopped by user before: %s", script))
      aborted_by_user <- TRUE
      break
    }
  }

  log_file <- file.path(
    "output/logs",
    sprintf("%s_%s_%s.log", timestamp, tools::file_path_sans_ext(basename(script)), mode)
  )
  tmp_log <- tempfile(pattern = "scspill_log_", fileext = ".txt")

  status <- system2(
    "Rscript",
    c(script),
    env = c(paste0("SCSPILL_MODE=", mode)),
    stdout = tmp_log,
    stderr = tmp_log
  )

  header <- build_log_header(
    script = script,
    mode = mode,
    target = target,
    lock_df = lock_df,
    seed_contract = script_seed_contract[[script]]
  )
  body <- if (file.exists(tmp_log)) readLines(tmp_log, warn = FALSE) else character(0)
  writeLines(c(header, body), con = log_file)

  if (status != 0) {
    stop(sprintf("Script failed: %s (see %s)", script, log_file))
  }
  message(sprintf("Completed: %s", script))
  message(sprintf("Log: %s", log_file))
}

if (!aborted_by_user) {
  check_expected_outputs(target)

  if (mode == "smoke") {
    message("Smoke run completed: pipeline integrity verified.")
  } else {
    message("Full run completed: manuscript-scale outputs generated.")
  }
  message("Replication workflow completed.")
} else {
  message("Partial interactive run completed. Expected-output checks were skipped.")
}
