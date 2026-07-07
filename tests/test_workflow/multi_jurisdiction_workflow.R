## Multi-jurisdiction (hierarchical/partially-pooled) workflow
##
## A fully self-contained example fitting two jurisdictions jointly, using
## fabricated data. Demonstrates:
##   - jurisdiction is only ever a dimension at the stack_jurisdictions()
##     step: each jurisdiction's data is prepared independently via
##     define_observation_data()/define_observation_model() (exactly like the
##     single-jurisdiction workflow), then combined via stack_jurisdictions(),
##     named by jurisdiction, before being passed to fit_waves()
##   - jurisdictions can have different/staggered date coverage -- every
##     per-jurisdiction vector is aligned to the same shared
##     target_infection_dates axis, so this is safe
##   - partial pooling: create_infection_timeseries() shares GP kernel
##     hyperparameters across jurisdictions (independent draws, shared
##     lengthscale/variance), and create_dow_priors() shares a hierarchical
##     day-of-week prior (nation-wide effect + per-jurisdiction deviation)
##     when dow_model = TRUE
##
## For a single-jurisdiction example, see single_jurisdiction_workflow.R.

devtools::load_all()
library(epiwave.params)
library(distributional)

set.seed(7)

target_infection_dates <- seq(as.Date("2024-01-01"), by = "day", length.out = 150)
cases_delay <- as_discrete_pmf(distributional::dist_gamma(shape = 3, rate = 0.5))

make_jurisdiction_cases <- function(start, end, lambda, dow_model) {
  dates <- target_infection_dates[start:end]
  counts <- rpois(length(dates), lambda = lambda)
  define_observation_data(
    timeseries_data = data.frame(date = dates, value = counts),
    delay_from_infection = cases_delay,
    proportion_infections = 0.5,
    dow_model = dow_model)
}

# Jurisdiction A: cases observed days 10-90.
# Jurisdiction B: cases observed days 60-140.
# Deliberately staggered/non-identical coverage, to exercise the
# master-date-axis alignment that makes joint fitting safe.
observation_model_A <- define_observation_model(
  target_infection_dates = target_infection_dates,
  cases = make_jurisdiction_cases(10, 90, lambda = 50, dow_model = TRUE)
)
observation_model_B <- define_observation_model(
  target_infection_dates = target_infection_dates,
  cases = make_jurisdiction_cases(60, 140, lambda = 30, dow_model = TRUE)
)

# named arguments to stack_jurisdictions(), named by jurisdiction -- this is
# the only place jurisdiction becomes a dimension
stacked <- stack_jurisdictions(
  jurisdiction_a = observation_model_A,
  jurisdiction_b = observation_model_B
)

fit_object <- fit_waves(
  observations = stacked,
  infection_model_type = "gp_growth_rate",
  n_chains = 2,
  warmup = 200,
  n_samples = 200
)

stopifnot(identical(dim(fit_object$infection_model), c(150L, 2L)))
stopifnot(identical(fit_object$jurisdictions, c("jurisdiction_a", "jurisdiction_b")))

## -- sanity check: staggered coverage stays correctly aligned ---------------
##
## A date only jurisdiction A has data for should show up as non-NA for A and
## NA for B in the underlying case matrix, and vice versa -- confirming row i
## means the same calendar date for every jurisdiction's column.

row_only_a <- 20   # within A's 10-90 window, outside B's 60-140 window
row_only_b <- 135  # within B's window, outside A's

stopifnot(!is.na(observation_model_A$observation_model_data$cases$case_vec[row_only_a]))
stopifnot(is.na(observation_model_B$observation_model_data$cases$case_vec[row_only_a]))
stopifnot(is.na(observation_model_A$observation_model_data$cases$case_vec[row_only_b]))
stopifnot(!is.na(observation_model_B$observation_model_data$cases$case_vec[row_only_b]))

rhats <- coda::gelman.diag(fit_object$fit, autoburnin = FALSE, multivariate = FALSE)
cat("max rhat:", max(rhats$psrf[, 1], na.rm = TRUE), "\n")

cat("\nMulti-jurisdiction workflow complete -- jurisdictions:",
    fit_object$jurisdictions, "\n")
