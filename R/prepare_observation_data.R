#' Prepare observation data
#'
#' @description Prepare the data objects needed for one jurisdiction's
#'  observation model for a single data stream (e.g. cases), called
#'  internally by `define_observation_model()` once per stream. Multiple
#'  jurisdictions are combined later, by `stack_jurisdictions()`.
#'
#' @param observation_data data for one jurisdiction, one stream, as
#'  returned by `define_observation_data()`.
#'  `timeseries_data` may be a plain data.frame/tibble with `date` and
#'  `value` columns -- it doesn't need to be pre-classed as
#'  `epiwave_timeseries`, see `as_epiwave_timeseries()`.
#'  `delay_from_infection` may be a single `discrete_pmf` object (replicated
#'  across `target_infection_dates`), or an already time-varying
#'  `discrete_pmf_series` (aligned to `target_infection_dates`)
#' @param target_infection_dates full date sequence of infection timeseries,
#'  shared across all jurisdictions in a fit
#'
#' @return a list describing this jurisdiction's stream: `convolution_matrix`,
#'  `case_vec`, `prop_vec`, `delays`, `dow_model`. Does not include GAM-based
#'  initial values -- those are only needed for `infection_model_type =
#'  'flat_prior'`, so they're computed lazily by `compute_flat_prior_inits()`
#'  when `fit_waves()` actually needs them, rather than unconditionally here.
#' @noRd
prepare_observation_data <- function (observation_data,
                                      target_infection_dates) {

  delays <- observation_data$delay_from_infection
  if (inherits(delays, 'discrete_pmf_series')) {
    delays <- delays[as.Date(target_infection_dates)]
    if (!identical(as.Date(delays$index), as.Date(target_infection_dates))) {
      stop('`delay_from_infection` dates must match `target_infection_dates`')
    }
  } else if (inherits(delays, 'discrete_pmf')) {
    delays <- epiwave.params::new_discrete_series(
      values = delays,
      index = target_infection_dates)
  } else {
    stop('`delay_from_infection` must be a discrete_pmf or ',
         'discrete_pmf_series object')
  }

  # a bare numeric prop (scalar or a vector already in date order) needs no
  # coercion at all -- as_matrix.numeric() recycles/validates it directly.
  # Only already-classed epiwave_timeseries objects (including
  # epiwave_greta_timeseries) need their dates checked against
  # target_infection_dates first.
  prop <- observation_data$proportion_infections
  if (inherits(prop, 'epiwave_timeseries')) {
    # epiwave_greta_timeseries stores dates nested under $timeseries$date
    # (it wraps a greta array alongside the date tibble), rather than a
    # flat $date column like epiwave_timeseries/epiwave_fixed_timeseries
    prop_dates <- if (inherits(prop, 'epiwave_greta_timeseries')) {
      prop$timeseries$date
    } else {
      prop$date
    }
    if (!identical(as.Date(prop_dates), as.Date(target_infection_dates))) {
      stop('`proportion_infections` dates must match `target_infection_dates`')
    }
  }

  timeseries_data <- as_epiwave_timeseries(observation_data$timeseries_data)
  case_vec <- as_matrix(timeseries_data, target_infection_dates)
  prop_vec <- as_matrix(prop, target_infection_dates)

  convolution_matrix <- new_convolution_matrix(delays)

  list(convolution_matrix = convolution_matrix,
       case_vec = case_vec,
       prop_vec = prop_vec,
       delays = delays,
       dow_model = observation_data$dow_model)
}
