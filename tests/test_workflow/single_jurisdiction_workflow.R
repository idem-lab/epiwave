## Single-jurisdiction workflow
##
## A minimal, fully self-contained example of the epiwave workflow for one
## jurisdiction, using fabricated data (no external/not-synced files needed).
## Demonstrates the current API end to end:
##   - raw data.frames passed straight to define_observation_data() (no
##     manual class(x) <- c(...) needed -- see as_epiwave_timeseries())
##   - delay distributions built with epiwave.params (discrete_pmf, combined
##     via `+`/add_discrete()), including a time-varying discrete_pmf_series
##   - proportion_infections as both a fixed value and a greta-array-derived
##     value (the IHR-from-CHR pattern)
##   - define_observation_data() -> define_observation_model() -> fit_waves()
##
## For a multi-jurisdiction (hierarchical/partially-pooled) example, see
## multi_jurisdiction_workflow.R.

devtools::load_all()
library(epiwave.params)
library(distributional)

set.seed(42)

## -- date range -----------------------------------------------------------

target_infection_dates <- seq(as.Date("2024-01-01"), by = "day", length.out = 150)

## -- delay distributions ----------------------------------------------------

# incubation period + onset-to-notification delay, combined by convolution
# (the `+` operator on discrete_pmf objects convolves them -- equivalent to
# add_discrete(incubation, onset_to_notification))
incubation <- as_discrete_pmf(distributional::dist_weibull(shape = 1.83, scale = 4.93))
onset_to_notification <- as_discrete_pmf(distributional::dist_gamma(shape = 3, rate = 0.5))
cases_delay <- incubation + onset_to_notification

# hospitalisation delay: a single distributional object works directly too
hosp_delay <- as_discrete_pmf(distributional::dist_weibull(shape = 2.51, scale = 10.17))

## -- fabricated observation data --------------------------------------------

# case counts: a plain date/value data.frame -- no epiwave_timeseries class
# needed, define_observation_data() coerces it automatically
case_dates <- target_infection_dates[20:130]
case_counts <- rpois(length(case_dates), lambda = 50 + 20 * sin(seq_along(case_dates) / 10))
notif_dat <- data.frame(date = case_dates, value = case_counts)

# hospitalisation counts, a subset of infections some days later
hosp_dates <- target_infection_dates[30:120]
hosp_counts <- rpois(length(hosp_dates), lambda = 5)
hosp_dat <- data.frame(date = hosp_dates, value = hosp_counts)

## -- proportions -------------------------------------------------------------

# case ascertainment rate: a fixed proportion
car <- 0.5

# hospitalisation proportion: derived from CHR (case-hospitalisation rate)
# and CAR via a greta array, so it can be estimated jointly with the rest of
# the model rather than fixed
chr <- greta::uniform(0, 1)
ihr <- create_epiwave_greta_timeseries(
  dates = target_infection_dates,
  car = car,
  chr_prior = chr)

## -- observation model, one jurisdiction -------------------------------------

observation_model <- define_observation_model(
  target_infection_dates = target_infection_dates,

  cases = define_observation_data(
    timeseries_data = notif_dat,
    delay_from_infection = cases_delay,
    proportion_infections = car,
    dow_model = TRUE),

  hospitalisations = define_observation_data(
    timeseries_data = hosp_dat,
    delay_from_infection = hosp_delay,
    proportion_infections = ihr)
)

## -- fit ----------------------------------------------------------------

# a single jurisdiction's define_observation_model() output goes straight to
# fit_waves() -- no list()/naming step needed. For multiple jurisdictions,
# see multi_jurisdiction_workflow.R's explicit stack_jurisdictions() step.
#
# flat_prior is the fastest infection model, good for a quick check like this
# one; for a real analysis, 'gp_infections'/'gp_growth_rate'/
# 'gp_growth_rate_deriv' model infections as a Gaussian Process instead (see
# create_infection_timeseries()'s docs for how the three GP formulations
# differ). Warmup/n_samples are kept small here for speed, not convergence.
fit_object <- fit_waves(
  observations = observation_model,
  infection_model_type = "flat_prior",
  n_chains = 2,
  warmup = 200,
  n_samples = 200
)

stopifnot(identical(dim(fit_object$infection_model), c(150L, 1L)))

rhats <- coda::gelman.diag(fit_object$fit, autoburnin = FALSE, multivariate = FALSE)
cat("max rhat:", max(rhats$psrf[, 1], na.rm = TRUE), "\n")

infection_draws <- greta::calculate(
  fit_object$infection_model,
  values = fit_object$fit,
  nsim = 200)[[1]]
cat("posterior median infections range:",
    round(range(apply(infection_draws, 2, median))), "\n")

## -- advanced: a time-varying delay distribution -----------------------------
##
## delay_from_infection can also be an already time-varying
## discrete_pmf_series, e.g. representing a change in testing/reporting
## delay partway through the study period.

notification_delay_early <- as_discrete_pmf(distributional::dist_gamma(shape = 3, rate = 0.5))
notification_delay_late <- as_discrete_pmf(distributional::dist_gamma(shape = 2, rate = 0.6))

time_varying_cases_delay <- new_discrete_series(
  values = c(
    rep(list(notification_delay_early), 75),
    rep(list(notification_delay_late), 75)
  ),
  index = target_infection_dates
)

observation_model_time_varying <- define_observation_model(
  target_infection_dates = target_infection_dates,
  cases = define_observation_data(
    timeseries_data = notif_dat,
    delay_from_infection = time_varying_cases_delay,
    proportion_infections = car,
    dow_model = TRUE)
)

cat("\nSingle-jurisdiction workflow complete.\n")
