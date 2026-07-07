#' Bundle a case ascertainment rate and case-hospitalisation-rate prior
#'
#' @description Builds a greta-backed `proportion_infections` input (e.g.
#'   the implied hospitalisation rate, IHR = CAR x CHR) without needing a
#'   date axis up front. The greta array this eventually becomes has to be
#'   dimensioned to `target_infection_dates`, but that axis is emergent --
#'   derived later, from the data, by `stack_jurisdictions()` -- so the
#'   actual `car * chr_prior` multiplication and dimensioning is deferred
#'   until then, once the axis is known.
#'
#' @param car case ascertainment rate; a fixed numeric value
#' @param chr_prior a greta array representing the prior distribution for
#'   the case hospitalisation rate
#'
#' @return a list with class `greta_proportion` containing `car` and
#'   `chr_prior`, unresolved until `stack_jurisdictions()`
#' @export
as_greta_timeseries <- function(car,
                                chr_prior) {

  out <- list(car = car, chr_prior = chr_prior)
  class(out) <- c("greta_proportion", class(out))
  out
}

#' Coerce a date/value table to an epiwave_fixed_timeseries object
#'
#' @description Observation data (case counts, hospitalisation counts) is
#'   naturally a table of `date`/`value` pairs with real per-date
#'   observations, rather than a
#'   single value replicated across dates -- unlike delay distributions or
#'   proportions, which usually apply uniformly. This function takes such a
#'   table (with genuinely partial date coverage, e.g. missing weekends, is
#'   fine -- gaps are handled downstream by `as_matrix()`) and classes it as
#'   an `epiwave_fixed_timeseries` object, so it doesn't need to be
#'   hand-classed with `class(x) <- c(...)`.
#'
#'   Already-classed `epiwave_timeseries` objects and plain numeric
#'   values/vectors are returned unchanged, so this is safe to call
#'   unconditionally on anything accepted as observation data.
#'
#' @param data an `epiwave_timeseries` object, a numeric value/vector, or a
#'   data.frame/tibble with `date` and `value` columns
#'
#' @return an `epiwave_timeseries` or numeric object, ready for `as_matrix()`
#' @noRd
as_epiwave_timeseries <- function(data) {

  if (inherits(data, "epiwave_timeseries") || is.numeric(data)) {
    return(data)
  }

  if (!all(c("date", "value") %in% names(data))) {
    stop("`data` must have `date` and `value` columns to be used as ",
         "observation data")
  }

  data <- tibble::as_tibble(data[c("date", "value")])
  data$date <- as.Date(data$date)

  class(data) <- c("epiwave_fixed_timeseries", "epiwave_timeseries", class(data))
  data
}
