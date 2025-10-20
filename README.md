**NOTE: This repository is under construction.**

# scspill

Bayesian synthetic control with spillovers (SAR). Provides `sc_spillover()` which
returns posterior draws of weights and rho, and identified treatment/spillover effects.

## Install

```r
install.packages(c("devtools","roxygen2","testthat","rmarkdown","progress","ggplot2","Matrix"))
devtools::document()
devtools::build()
devtools::check()
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