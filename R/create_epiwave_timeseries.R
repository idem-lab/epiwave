#' Create a greta-compatible timeseries object
#'
#' @description The epiwave model functions expect data in a long format,
#'   structured to have a value for every date. This function creates a
#'   greta-compatible timeseries from a case ascertainment rate and a
#'   hospitalisation rate prior.
#'
#' @param dates infection dates sequence
#' @param car case ascertainment rate; either a numeric value or an
#'   `epiwave_timeseries` object
#' @param chr_prior a greta array representing the prior distribution for
#'   the case hospitalisation rate
#'
#' @importFrom tibble tibble
#'
#' @return a list with class `epiwave_greta_timeseries` containing
#'   components `timeseries` (a tibble of dates) and `ihr` (a greta array of
#'   implied hospitalisation rates)
#' @export
create_epiwave_greta_timeseries <- function(dates,
                                            car,
                                            chr_prior) {
  long_unique <- tibble::tibble(date = dates)

  dim(chr_prior) <- length(dates)

  if ("epiwave_timeseries" %in% class(car)) {
    car <- car$value
  }

  ihr_greta <- car * chr_prior

  long_combined <- list(timeseries = long_unique,
                        ihr = ihr_greta)

  class(long_combined) <- c("epiwave_greta_timeseries",
                            "epiwave_timeseries",
                            class(long_combined))
  long_combined
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
#' @export
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
