#' Prepare observation data
#'
#' @description Clean and validate the data needed for one jurisdiction's
#'  observation model for a single data stream (e.g. cases), called
#'  internally by `define_observation_model()` once per stream. Deliberately
#'  does *not* take a `target_infection_dates` argument: this layer keeps
#'  everything in its own native, date-keyed form (real `Date`s, not
#'  positions), and doesn't decide a shared date axis at all -- that's
#'  derived later, from the data, by `stack_jurisdictions()`. This is what
#'  lets `target_infection_dates` be emergent rather than something the user
#'  has to compute correctly up front.
#'
#' @param observation_data data for one jurisdiction, one stream, as
#'  returned by `define_observation_data()`.
#'  `timeseries_data` may be a plain data.frame/tibble with `date` and
#'  `value` columns -- it doesn't need to be pre-classed, see
#'  `as_epiwave_timeseries()`.
#'  `delay_from_infection` may be a single `discrete_pmf` object, or an
#'  already time-varying `discrete_pmf_series`
#'
#' @return a list describing this jurisdiction's stream: `timeseries_data`
#'  (a date-keyed `epiwave_timeseries` object), `delays` (native
#'  `discrete_pmf`/`discrete_pmf_series`), `proportion_infections` (native,
#'  unresolved if greta-backed), `dow_model`, and `implied_range` -- a
#'  length-2 `Date` vector giving the earliest and latest infection dates
#'  this stream's data could inform (earliest = first observed date minus
#'  the delay distribution's maximum support; latest = last observed date).
#'  Does not build a convolution matrix or align anything to a shared axis
#'  yet -- both need `target_infection_dates`, which isn't decided until
#'  `stack_jurisdictions()`.
#' @noRd
prepare_observation_data <- function (observation_data) {

  delays <- observation_data$delay_from_infection
  if (!inherits(delays, c('discrete_pmf', 'discrete_pmf_series'))) {
    stop('`delay_from_infection` must be a discrete_pmf or ',
         'discrete_pmf_series object')
  }

  timeseries_data <- as_epiwave_timeseries(observation_data$timeseries_data)
  if (!inherits(timeseries_data, 'epiwave_timeseries')) {
    stop('`timeseries_data` must resolve to a date/value table -- ',
         'a data.frame/tibble with `date` and `value` columns, or an ',
         'already-classed epiwave_timeseries object')
  }
  observed_dates <- as.Date(timeseries_data$date)

  # proportion_infections' *date-alignment* can't be checked yet -- that
  # needs target_infection_dates, which isn't decided until
  # stack_jurisdictions() (see align_stream_to_axis()) -- but a basic type
  # check doesn't need the axis, so it can and should still happen here,
  # rather than only surfacing once fit_waves()/stack_jurisdictions() runs
  prop <- observation_data$proportion_infections
  if (!is.numeric(prop) &&
      !inherits(prop, c('epiwave_timeseries', 'greta_proportion'))) {
    stop('`proportion_infections` must be numeric, an epiwave_timeseries ',
         'object, or a greta_proportion object (see as_greta_timeseries())')
  }

  max_delay <- if (inherits(delays, 'discrete_pmf_series')) {
    max(vapply(delays$values, function(x) max(x$step), numeric(1)))
  } else {
    max(delays$step)
  }

  implied_range <- c(
    min(observed_dates) - max_delay,
    max(observed_dates)
  )

  # a discrete_pmf_series was built by the caller with its own fixed index,
  # before target_infection_dates existed to build it against -- if that
  # index doesn't even cover this stream's own implied range, stacking will
  # fail later regardless of what any other stream/jurisdiction needs, so
  # fail now with a message that says what's wrong and how to fix it,
  # rather than surfacing as a cryptic subsetting error deep in fit_waves().
  # (Combining with other streams/jurisdictions can *still* widen the final
  # axis beyond what's checked here -- that can only be caught once
  # everything is known, in stack_jurisdictions()/align_stream_to_axis().)
  if (inherits(delays, 'discrete_pmf_series')) {
    series_range <- range(as.Date(delays$index))
    if (implied_range[1] < series_range[1] || implied_range[2] > series_range[2]) {
      stop(sprintf(
        paste0(
          "`delay_from_infection`'s discrete_pmf_series only covers %s to %s, ",
          "but this stream's own data implies infections need to be tracked ",
          "over %s to %s (the earliest/latest observed date, adjusted for the ",
          "delay distribution's reach). Extend the series' index to cover at ",
          "least this range -- note the final fitting axis, once combined with ",
          "any other streams/jurisdictions, may need to extend further still."),
        format(series_range[1]), format(series_range[2]),
        format(implied_range[1]), format(implied_range[2])))
    }
  }

  list(
    timeseries_data = timeseries_data,
    delays = delays,
    proportion_infections = observation_data$proportion_infections,
    dow_model = observation_data$dow_model,
    implied_range = implied_range
  )
}
