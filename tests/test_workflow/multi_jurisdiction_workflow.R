## Multi-jurisdiction (hierarchical/partially-pooled) workflow
##
## A fully self-contained example fitting two jurisdictions jointly, using
## fabricated data. Demonstrates:
##   - jurisdiction is only ever a dimension at the stack_jurisdictions()
##     step: each jurisdiction's data is prepared independently via
##     define_observation_data()/define_observation_model() (exactly like the
##     single-jurisdiction workflow, and just as emergent -- no
##     target_infection_dates supplied anywhere), then combined via
##     stack_jurisdictions(), named by jurisdiction, before being passed to
##     fit_waves()
##   - the combined date axis is the emergent union of every jurisdiction's
##     own implied range; jurisdictions can have different/staggered date
##     coverage, and every stream is aligned to the shared axis by matching
##     real dates (not by positional assumption), so this is safe
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

# a date range purely for fabricating synthetic data, not target_infection_dates
dates <- seq(as.Date("2024-01-01"), by = "day", length.out = 150)
cases_delay <- as_discrete_pmf(distributional::dist_gamma(shape = 3, rate = 0.5))

# Jurisdiction A: cases observed days 10-90.
# Jurisdiction B: cases observed days 60-140.
# Deliberately staggered/non-identical coverage, to exercise the
# master-date-axis alignment that makes joint fitting safe.
case_dates_a <- dates[10:90]
case_counts_a <- rpois(length(case_dates_a), lambda = 50)
observation_model_A <- define_observation_model(
  cases = define_observation_data(
    timeseries_data = data.frame(date = case_dates_a, value = case_counts_a),
    delay_from_infection = cases_delay,
    proportion_infections = 0.5,
    dow_model = TRUE)
)

case_dates_b <- dates[60:140]
case_counts_b <- rpois(length(case_dates_b), lambda = 30)
observation_model_B <- define_observation_model(
  cases = define_observation_data(
    timeseries_data = data.frame(date = case_dates_b, value = case_counts_b),
    delay_from_infection = cases_delay,
    proportion_infections = 0.5,
    dow_model = TRUE)
)

# named arguments to stack_jurisdictions(), named by jurisdiction -- this is
# the only place jurisdiction becomes a dimension, and where the emergent
# target_infection_dates axis (the union of both jurisdictions' implied
# ranges) is actually derived
stacked <- stack_jurisdictions(
  jurisdiction_a = observation_model_A,
  jurisdiction_b = observation_model_B
)

# worth inspecting the derived axis and alignment before committing to a
# (potentially slow) fit -- print() for a quick text summary,
# plot_observation_coverage() to see it
print(stacked)
plot_observation_coverage(stacked)

fit_object <- fit_waves(
  observations = stacked,
  infection_model_type = "gp_growth_rate",
  n_chains = 2,
  warmup = 200,
  n_samples = 200
)

stopifnot(ncol(fit_object$infection_model) == 2L)
stopifnot(identical(fit_object$jurisdictions, c("jurisdiction_a", "jurisdiction_b")))

cat("emergent fitting window:", format(range(fit_object$infection_days)),
    "(", length(fit_object$infection_days), "days )\n")

## -- sanity check: staggered coverage stays correctly aligned ---------------
##
## A date only jurisdiction A has data for should show up as non-NA for A and
## NA for B in the stacked case matrix, and vice versa -- confirming row i
## means the same calendar date for every jurisdiction's column, found by
## matching real dates rather than assuming fixed row numbers (the emergent
## axis's start depends on the delay distribution's reach, not just the
## fabrication window above).

case_mat <- stacked$observation_model_data$cases$case_mat
row_only_a <- match(dates[20], stacked$target_infection_dates)   # within A's 10-90 window, outside B's 60-140 window
row_only_b <- match(dates[135], stacked$target_infection_dates)  # within B's window, outside A's

stopifnot(!is.na(case_mat[row_only_a, "jurisdiction_a"]))
stopifnot(is.na(case_mat[row_only_a, "jurisdiction_b"]))
stopifnot(is.na(case_mat[row_only_b, "jurisdiction_a"]))
stopifnot(!is.na(case_mat[row_only_b, "jurisdiction_b"]))

rhats <- coda::gelman.diag(fit_object$fit, autoburnin = FALSE, multivariate = FALSE)
cat("max rhat:", max(rhats$psrf[, 1], na.rm = TRUE), "\n")

cat("\nMulti-jurisdiction workflow complete -- jurisdictions:",
    fit_object$jurisdictions, "\n")
