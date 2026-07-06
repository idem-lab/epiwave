
# epiwave

This package is the first implementation of a novel algorithm to
estimate infection incidence timeseries. It relies underneath on
existing sampling algorithms from the `greta` software. we get infection
incidence, no other package does this. estimating r effective. we could
benchmark but pred accuracy not super important for users.

## Package to estimate epidemiological parameters using a Gaussian Process for infection timeseries

[![source](https://img.shields.io/badge/source-GitHub-success?style=flat&labelColor=gray)](https://github.com/idem-lab/epiwave)
[![wip](https://www.repostatus.org/badges/latest/wip.svg)](https://www.repostatus.org/#wip)
[![v.0.1.0-alpha](https://zenodo.org/badge/720975798.svg)](https://zenodo.org/doi/10.5281/zenodo.10398079)

`epiwave` estimates a time-series of the daily number of new infections,
using a time-series of daily reported case numbers and an estimated
distribution of the time delay between the onset of infection and case
notification. This model is unique because it integrates multiple data
sources to estimate the underlying infection dynamics. Besides case,
hospitalisation, and death notification data, these data sources could
include wastewater data and serological data. From the estimated
time-series of daily new infections we can then estimate useful epidemic
parameters such as the effective reproduction number over time. This
package will be especially useful when case incidence data, the most
commonly used source of information for estimating effective
reproduction number, is unreliable or challenging to use, but other data
sources are available.

## Installation

`epiwave` can be installed from [GitHub](https://github.com/) with:

``` r
remotes::install_github('idem-lab/epiwave')
```

## Usage

### One jurisdiction

``` r
library(epiwave)
library(epiwave.params)
library(distributional)

dates <- seq(as.Date("2024-01-01"), by = "day", length.out = 150)

# delay from infection to case notification
cases_delay <- as_discrete_pmf(distributional::dist_gamma(shape = 3, rate = 0.5))

# case counts: a plain date/value data.frame, no need to pre-class it
notif_dat <- data.frame(
  date = dates[20:130],
  value = rpois(111, lambda = 50)
)

observation_model <- define_observation_model(
  target_infection_dates = dates,
  cases = define_observation_data(
    timeseries_data = notif_dat,
    delay_from_infection = cases_delay,
    proportion_infections = 0.5,
    dow_model = TRUE)
)

# a single jurisdiction goes straight to fit_waves(), no combining step needed
fit_object <- fit_waves(
  observations = observation_model,
  infection_model_type = "flat_prior"
)
```

### Multiple jurisdictions

Jurisdictions are prepared independently, then combined explicitly via
`stack_jurisdictions()`, named by jurisdiction. Jurisdictions sharing
one fit are partially pooled (shared GP kernel hyperparameters, and a
shared hierarchical day-of-week prior where requested).

``` r
observation_model_b <- define_observation_model(
  target_infection_dates = dates,
  cases = define_observation_data(
    timeseries_data = data.frame(date = dates[60:140], value = rpois(81, lambda = 30)),
    delay_from_infection = cases_delay,
    proportion_infections = 0.5,
    dow_model = TRUE)
)

stacked <- stack_jurisdictions(
  jurisdiction_a = observation_model,
  jurisdiction_b = observation_model_b
)

fit_object <- fit_waves(
  observations = stacked,
  infection_model_type = "gp_growth_rate"
)
```

## Citation

When using this package, please cite the underlying statistical
software, `greta` as well as the package itself:

## Contribution

This is a work in progress. If you see any mistakes in the package
(branch `main`), let us know by logging a bug on the
[issues](https://github.com/idem-lab/epiwave/issues) page.

## Code of Conduct

Please note that the `epiwave` project is released with a [Contributor
Code of
Conduct](https://idem-lab.github.io/epiwave/CODE_OF_CONDUCT.html). By
contributing to this project, you agree to abide by its terms.

## Support

This project was supported by ???? BDSS??? \[the Australia-Aotearoa
Consortium for Epidemic Forecasting and Analytics.\]

<a href="https://acefa-hubs.github.io/"><img src="man/figures/ACEFA.png" align = "center" height="150" alt="EpiStrainDynamics website" /></a>
