#' Prepare observation data
#'
#' @description Prepare the data objects needed for one jurisdiction's
#'  observation model for a single data stream (e.g. cases). Multiple
#'  jurisdictions are combined later, by `stack_jurisdictions()`.
#'
#' @param observation_data data for one jurisdiction, one stream, as
#'  returned by `define_observation_data()`/`define_sero_data()`.
#'  `delay_from_infection` may be a single `discrete_pmf`/`discrete_weights`
#'  object (replicated across `target_infection_dates`), or an already
#'  time-varying `discrete_pmf_series`/`discrete_weights_series` (aligned to
#'  `target_infection_dates`)
#' @param target_infection_dates full date sequence of infection timeseries,
#'  shared across all jurisdictions in a fit
#'
#' @return a list describing this jurisdiction's stream: `convolution_matrix`,
#'  `case_vec`, `prop_vec`, `dow_model`, `inits_values`, `observable_idx`, and
#'  (for seroprevalence streams) `total_pop`/`size_vec`
#' @export
prepare_observation_data <- function (observation_data,
                                      target_infection_dates) {

  delays <- observation_data$delay_from_infection
  if (inherits(delays, c('discrete_pmf_series', 'discrete_weights_series'))) {
    delays <- delays[as.Date(target_infection_dates)]
    if (!identical(as.Date(delays$index), as.Date(target_infection_dates))) {
      stop('`delay_from_infection` dates must match `target_infection_dates`')
    }
  } else if (inherits(delays, c('discrete_pmf', 'discrete_weights'))) {
    delays <- epiwave.params::new_discrete_series(
      values = delays,
      index = target_infection_dates)
  } else {
    stop('`delay_from_infection` must be a discrete_pmf, discrete_weights, ',
         'discrete_pmf_series, or discrete_weights_series object')
  }

  prop <- observation_data$proportion_infections
  if (!inherits(prop, 'epiwave_timeseries')) {
    prop <- create_epiwave_timeseries(
      dates = target_infection_dates,
      value = prop)
  } else {
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

  case_vec <- as_matrix(observation_data$timeseries_data, target_infection_dates)
  prop_vec <- as_matrix(prop, target_infection_dates)

  convolution_matrix <- new_convolution_matrix(delays)

  inits <- inits_by_jurisdiction(
    case_vec,
    delays,
    prop_vec,
    target_infection_dates)

  inits_values <- rep(NA_real_, length(target_infection_dates))
  observable_idx <- rep(FALSE, length(target_infection_dates))
  inits_values[inits$observable_idx] <- inits$inits_values
  observable_idx[inits$observable_idx] <- TRUE

  out <- list(convolution_matrix = convolution_matrix,
              case_vec = case_vec,
              prop_vec = prop_vec,
              dow_model = observation_data$dow_model,
              inits_values = inits_values,
              observable_idx = observable_idx)

  if ('total_pop' %in% names(observation_data)) {
    out$total_pop <- observation_data$total_pop
  }
  if ('size_vec' %in% names(observation_data)) {
    out$size_vec <- as_matrix(observation_data$size_vec, target_infection_dates)
  }

  out
}
