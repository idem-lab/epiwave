# epiwave
## Package to estimate epidemiological parameters using a Gaussian Process for infection timeseries

[![source](https://img.shields.io/badge/source-GitHub-success?style=flat&labelColor=gray)](https://github.com/idem-lab/epiwave)
[![wip](https://www.repostatus.org/badges/latest/wip.svg)](https://www.repostatus.org/#wip)
[![v.0.1.0-alpha](https://zenodo.org/badge/720975798.svg)](https://zenodo.org/doi/10.5281/zenodo.10398079)

`epiwave` estimates a time-series of the daily number of new infections, using a time-series of daily reported case numbers and an estimated distribution of the time delay between the onset of infection and case notification.
This model is unique because it integrates multiple data sources to estimate the underlying infection dynamics.
Besides case, hospitalisation, and death notification data, these data sources could include wastewater data and serological data.
From the estimated time-series of daily new infections we can then estimate useful epidemic parameters such as the effective reproduction number over time.
This package will be especially useful when case incidence data, the most commonly used source of information for estimating effective reproduction number, is unreliable or challenging to use, but other data sources are available.

`epiwave` can be installed with GitHub, as follows:
``` r
remotes::install_github('idem-lab/epiwave')
```
