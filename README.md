**NOTE: This repository is under construction.**

# scspill

Bayesian synthetic control with spillovers (SAR). Provides `sc_spillover()` which
returns posterior draws of weights and rho, and identified treatment/spillover effects.

## Quick Start Guide

This replication package contains all code and data necessary to reproduce the analyses presented in the paper.

### System Requirements
- **R version**: >= 4.2.0
- **Operating system**: Windows, macOS, or Linux
- **C++ compiler**: Required for building C++ components (automatically handled by R on most systems)

### Quick Reproduction Steps

1. **Install required R packages:**
```r
# Core dependencies
install.packages(c("Matrix", "progress", "ggplot2", "Rcpp", "RcppArmadillo"))

# Additional packages for vignettes
install.packages(c("dplyr", "tidyr", "knitr", "rmarkdown", "kableExtra", "ggpattern"))

# Development tools (if not already installed)
install.packages(c("devtools", "roxygen2"))
```

2. **Load the package source code:**
```r
# Set working directory to the package root
setwd("/path/to/scspillover")

# Load all functions from source
devtools::load_all(".")
```

3. **Run the vignettes to reproduce results:**
```r
# Option A: Using RStudio
# Open each .Rmd file in vignettes/ and click "Knit"

# Option B: Using command line
rmarkdown::render("vignettes/vignette-california.Rmd")
rmarkdown::render("vignettes/vignette-sudan.Rmd")
rmarkdown::render("vignettes/vignette-simulation.Rmd")
rmarkdown::render("vignettes/vignette-geweke.Rmd")
```

### What Each Vignette Reproduces

- **`vignette-california.Rmd`**: Main application - California tobacco tax analysis (Tables and Figures)
- **`vignette-sudan.Rmd`**: Main application - Sudan split analysis (Tables and Figures)
- **`vignette-simulation.Rmd`**: Monte Carlo simulation study results
- **`vignette-geweke.Rmd`**: MCMC convergence diagnostics

**Note**: All data files needed for reproduction are included in the `data/` directory. No external data downloads are required.

### Expected Runtime

- California vignette: ~30-60 minutes (depending on hardware)
- Sudan vignette: ~30-60 minutes
- Simulation vignette: ~2-4 hours (multiple scenarios)
- Geweke vignette: ~10-20 minutes

### Troubleshooting

- If C++ compilation fails, ensure you have a C++ compiler installed (Rtools on Windows, Xcode Command Line Tools on macOS, build-essential on Linux)
- If vignettes fail to render, ensure all required packages listed above are installed
- For issues with specific vignettes, check that the working directory is set to the package root

## Required Packages

### Core Dependencies
- `stats`, `methods`: Base R packages
- `Matrix`: Sparse matrix operations
- `progress`: Progress bars
- `ggplot2`: Plotting
- `Rcpp` (>= 1.0.12): R-C++ integration
- `RcppArmadillo`: Linear algebra library (C++)

### Additional Packages for Vignettes
- `dplyr`, `tidyr`: Data manipulation
- `knitr`, `rmarkdown`: RMarkdown document processing
- `kableExtra`: Table formatting
- `ggpattern`: Graph patterns

### Development/Build Tools
- `devtools`, `roxygen2`: Package development
- `testthat` (>= 3.1.0): Testing framework
- `covr`, `pkgdown`: Code coverage and documentation site generation

## Directory Structure

### `data/`
Contains data files used for analysis.
- `raw/`: Raw data files (from World Bank, IMF, etc.)
- `processed/`: Preprocessed data (panel data, weight matrices, etc.)
- `.rda` files: R data objects (California tobacco, Sudan split datasets, etc.)

### `R/`
R source files implementing the package's main functionality.
- `10_sc_spillover.R`: Main function `sc_spillover()` implementation
- `21_mcmc_alpha.R`, `22_mcmc_sar.R`: MCMC sampling routines
- `01_utils.R` - `04_utils_*.R`: Utility functions (data preparation, plotting, diagnostics)
- `30_simulation_.R`: Simulation study functions
- `40_geweke_latest.R`, `41_robustness_check.R`: Diagnostics and robustness checks
- `61_processed_california.R`, `62_processed_sudan.R`: Data preprocessing scripts

### `src/`
C++ source files for performance-critical computations (using Rcpp).
- `20_mcmc.cpp`: MCMC sampler C++ implementation
- `40_geweke_latest.cpp`: Geweke diagnostics C++ implementation

### `vignettes/`
RMarkdown files for reproducing the paper's analyses.
- `vignette-california.Rmd`: California tobacco tax analysis
- `vignette-sudan.Rmd`: Sudan split analysis
- `vignette-simulation.Rmd`: Monte Carlo simulation study
- `vignette-geweke.Rmd`: Geweke diagnostics analysis

## Installation

```r
# Install development dependencies
install.packages(c("devtools", "roxygen2", "testthat", "rmarkdown"))

# Install package dependencies
install.packages(c("Matrix", "progress", "ggplot2", "Rcpp", "RcppArmadillo"))

# Build and install the package
devtools::clean_dll()
devtools::document()
devtools::build()
devtools::load_all(".")
devtools::check()
```

For local installation:

```r
tmp_lib <- file.path(tempdir(), "lib_scspill")
remotes::install_local(
  path = ".",
  lib = tmp_lib,
  build_vignettes = FALSE,
  upgrade = "never",
  INSTALL_opts = c("--no-multiarch"),
  force = TRUE
)
```

## Usage

### Basic Example

```r
library(scspill)
# panel: data.frame(unit, time, y)
# w: numeric vector (length N), W: row-normalized NxN matrix
fit <- sc_spillover(panel, treated_unit = "CA", T0 = 18, w = w, W = W,
                    M = 2000, burn = 1000, seed = 1)
summary(fit)
plot(fit, type = "treated")
```

### Running Vignettes

To reproduce the analyses from the paper, run the RMarkdown files in the `vignettes/` directory.

**Option 1: Using RStudio**
- Open each vignette file (`.Rmd`) in RStudio
- Click the "Knit" button

**Option 2: Using R command line**

```r
rmarkdown::render("vignettes/vignette-california.Rmd")
rmarkdown::render("vignettes/vignette-sudan.Rmd")
rmarkdown::render("vignettes/vignette-simulation.Rmd")
rmarkdown::render("vignettes/vignette-geweke.Rmd")
```

Available vignettes:
- `vignette-california.Rmd`: California tobacco tax application
- `vignette-sudan.Rmd`: Sudan split application
- `vignette-simulation.Rmd`: Monte Carlo simulation study
- `vignette-geweke.Rmd`: Geweke diagnostics

## Citation

If you use this package, please cite the accompanying paper: XXX

## Data Sources

- World Bank DataBank: https://databank.worldbank.org/source/world-development-indicators#
- IMF: https://data.imf.org/en/Data-Explorer?datasetUrn=IMF.STA:IMTS(1.0.0) or https://data.imf.org/en/datasets/IMF.STA:IMTS
