#' Expand constant value into a long tibble
#'
#' @description The epiwave model functions expect data in a long format,
#'   structured to have a value for every date. This function creates a
#'   tibble of this structure from a single value that is replicated in each
#'   cell.
#'
#' @param dates infection dates sequence
#' @param value value to be replicated in each cell
#'
#' @importFrom tibble tibble
#'
#' @return a long tibble with the constant value replicated across all
#'   dates
#' @export
create_epiwave_timeseries <- function(dates,
                                      value) {
  long_combined <- tibble::tibble(date = dates,
                                  value = value)

  class(long_combined) <- c("epiwave_timeseries", class(long_combined))
  long_combined
}

#' Expand a discrete PMF into a long tibble
#'
#' @description The epiwave model functions expect data in a long format,
#'   structured to have a value for every date. This function creates a
#'   tibble of this structure from a single `discrete_pmf` or
#'   `discrete_pmf_series` object that is replicated in each cell.
#'
#' @param dates infection dates sequence
#' @param value a `discrete_pmf` or `discrete_pmf_series` object to be
#'   replicated in each cell
#'
#' @return a long tibble with the distribution replicated across all
#'   dates, with class `epiwave_massfun_timeseries`
#' @export
create_epiwave_massfun_timeseries <- function(dates,
                                              value) {
  timeseries <- create_epiwave_timeseries(
    dates = dates,
    value = list(value)
  )

  class(timeseries) <- c("epiwave_massfun_timeseries",
                         class(value)[1],
                         class(timeseries))
  timeseries
}

#' Expand a fixed value into a long tibble
#'
#' @description The epiwave model functions expect data in a long format,
#'   structured to have a value for every date. This function creates a
#'   tibble of this structure from a single fixed numeric value that is
#'   replicated in each cell.
#'
#' @param dates infection dates sequence
#' @param value a fixed numeric value to be replicated in each cell
#'
#' @return a long tibble with the fixed value replicated across all
#'   dates, with class `epiwave_fixed_timeseries`
#' @export
create_epiwave_fixed_timeseries <- function(dates,
                                            value) {
  timeseries <- create_epiwave_timeseries(
    dates = dates,
    value = value
  )

  class(timeseries) <- c("epiwave_fixed_timeseries", class(timeseries))
  timeseries
}

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
