**NOTE: This repository is under construction.**

# scspill

Bayesian synthetic control with spillovers (SAR). Provides `sc_spillover()` which
returns posterior draws of weights and rho, and identified treatment/spillover effects.

## Install

```r
install.packages(c("devtools","roxygen2","testthat","rmarkdown","progress","ggplot2","Matrix"))
devtools::clean_dll()
devtools::document()
devtools::build()
devtools::load_all(".")
devtools::check()

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

```r
# remotes::install_github("yourname/scspill")
library(scspill)
# panel: data.frame(unit, time, y)
# w: numeric vector (length N), W: row-normalized NxN matrix
fit <- sc_spillover(panel, treated_unit = "CA", T0 = 18, w = w, W = W,
                    M = 2000, burn = 1000, seed = 1)
summary(fit)
plot(fit, type = "treated")
```

## Replication Package
- Vignettes reproduce California and Sudan applications.
- To rebuild docs: devtools::document(); to check: devtools::check().
- For a site: pkgdown::build_site().

## Citation

If you use this package, please cite the accompanying paper: XXX


Data
- World Bank DataBank: https://databank.worldbank.org/source/world-development-indicators#
- IMF: https://data.imf.org/en/Data-Explorer?datasetUrn=IMF.STA:IMTS(1.0.0) or https://data.imf.org/en/datasets/IMF.STA:IMTS
