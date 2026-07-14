#' Define observation model
#'
#' @description
#' Bundle one jurisdiction's observation streams (e.g. cases,
#' hospitalisations) together. Call this once per jurisdiction. A
#' single-jurisdiction fit can pass this straight to `fit_waves()`; to
#' combine multiple jurisdictions, pass several of these to
#' `stack_jurisdictions()` first.
#'
#' No `target_infection_dates` argument is needed here -- the date axis is
#' emergent, derived later (by `stack_jurisdictions()`) from every stream's
#' own data and delay distributions. This function just aggregates each
#' stream's implied date range into one range for the whole jurisdiction.
#'
#' @param ... observation data sets for this jurisdiction (as returned by
#'  `define_observation_data()`), named by stream
#'
#' @return list describing one jurisdiction's observation model, with class
#'  `epiwave_observation_model`
#' @export
#'
define_observation_model <- function (...) {

  observation_list <- list(...)

  prepared_observation_model_data <- lapply(observation_list,
                                            prepare_observation_data)

  implied_ranges <- lapply(prepared_observation_model_data,
                           function(x) x$implied_range)
  starts <- do.call(c, lapply(implied_ranges, `[`, 1))
  ends <- do.call(c, lapply(implied_ranges, `[`, 2))
  implied_range <- c(min(starts), max(ends))

  out <- list(observation_model_data = prepared_observation_model_data,
             implied_range = implied_range)

  class(out) <- c("epiwave_observation_model", class(out))
  out

}
