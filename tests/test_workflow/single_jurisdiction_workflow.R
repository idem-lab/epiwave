## Single-jurisdiction workflow
##
## A minimal, fully self-contained example of the epiwave workflow for one
## jurisdiction, using fabricated data (no external/not-synced files needed).
## Demonstrates the current API end to end:
##   - raw data.frames passed straight to define_observation_data() (no
##     manual class(x) <- c(...) needed -- no need to pre-class)
##   - delay distributions built with epiwave.params (discrete_pmf, combined
##     via `+`/add_discrete()), including a time-varying discrete_pmf_series
##   - proportion_infections as both a fixed value and a greta-array-derived
##     value (the IHR-from-CHR pattern)
##   - define_observation_data() -> define_observation_model() -> fit_waves()
##   - target_infection_dates is emergent: it's never supplied anywhere in
##     this script -- fit_waves() derives it from the data and the delay
##     distributions (earliest = first observed date minus the longest
##     delay's reach; latest = last observed date)
##
## For a multi-jurisdiction (hierarchical/partially-pooled) example, see
## multi_jurisdiction_workflow.R.

devtools::load_all()
library(epiwave.params)
library(distributional)

set.seed(42)

## -- a date range purely for fabricating synthetic data ---------------------
##
## This is NOT target_infection_dates -- it's just a convenient window to
## generate example data over. The actual fitting date axis is derived later,
## and will typically extend a bit earlier than any of this (to cover the
## delay distributions' reach).

dates <- seq(as.Date("2024-01-01"), by = "day", length.out = 150)

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
case_dates <- dates[20:130]
case_counts <- rpois(length(case_dates), lambda = 50 + 20 * sin(seq_along(case_dates) / 10))
notif_dat <- data.frame(date = case_dates, value = case_counts)

# hospitalisation counts, a subset of infections some days later
hosp_dates <- dates[30:120]
hosp_counts <- rpois(length(hosp_dates), lambda = 5)
hosp_dat <- data.frame(date = hosp_dates, value = hosp_counts)

## -- proportions -------------------------------------------------------------

# case ascertainment rate: a fixed proportion
car <- 0.5

# hospitalisation proportion: derived from CHR (case-hospitalisation rate)
# and CAR via a greta array, so it can be estimated jointly with the rest of
# the model rather than fixed. No dates needed here -- the greta array can
# only be dimensioned once the fitting axis is known, so that's deferred
# until fit_waves()/stack_jurisdictions() actually needs it.
chr <- greta::uniform(0, 1)
ihr <- as_greta_timeseries(car = car, chr_prior = chr)

## -- observation model, one jurisdiction -------------------------------------

observation_model <- define_observation_model(

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

stopifnot(ncol(fit_object$infection_model) == 1L)

cat("emergent fitting window:", format(range(fit_object$infection_days)),
    "(", length(fit_object$infection_days), "days )\n")
cat("earliest observed data: cases", format(min(case_dates)),
    ", hospitalisations", format(min(hosp_dates)), "\n")

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
## delay partway through the study period. Since the fitting axis is
## emergent, a time-varying series' own index needs to cover whatever that
## axis turns out to be -- build it over a generously wide range (here, the
## same window used to fabricate data above comfortably covers it).

notification_delay_early <- as_discrete_pmf(distributional::dist_gamma(shape = 3, rate = 0.5))
notification_delay_late <- as_discrete_pmf(distributional::dist_gamma(shape = 2, rate = 0.6))

time_varying_cases_delay <- new_discrete_series(
  values = c(
    rep(list(notification_delay_early), 75),
    rep(list(notification_delay_late), 75)
  ),
  index = dates
)

observation_model_time_varying <- define_observation_model(
  cases = define_observation_data(
    timeseries_data = notif_dat,
    delay_from_infection = time_varying_cases_delay,
    proportion_infections = car,
    dow_model = TRUE)
)

cat("\nSingle-jurisdiction workflow complete.\n")
