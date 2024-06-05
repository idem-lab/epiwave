#' Expand constant value into long tibble
#'
#' @description The lowerGPreff model functions expect data in a long format,
#'  which is structed to have a value for every unique date and jurisdiction
#'  pair. This function create a tibble of this structure out of a single
#'  value that should be replicated in each cell.
#'
#' @param dates infection dates sequence
#' @param jurisdictions jurisdiction names
#' @param value value to be replicated in each cell
#'
#' @importFrom dplyr mutate
#' @importFrom tibble tibble
#'
#' @return long tibble with constant value replicated
#' @export
#Assign superclass structure for timeseries data (i.e. distributions)
create_lowerGPreff_timeseries <- function (dates,
                                           jurisdictions,
                                           value) {

  long_unique <- expand.grid(date = dates,
                             jurisdiction = jurisdictions)
  long_combined <- long_unique |>
    dplyr::mutate(value = value) |> #consider another name for constant_val
    #dplyr::mutate(!!col_name := constant_val) |>
    tibble::tibble()

  class(long_combined) <- c("lowerGPreff_timeseries", class(long_combined))

  long_combined

  #also need another function that creates a proper time varying time series
  #(i.e. not outputting constant value)
  #:= name injection, assigns variable name to column name
}

#' Expand distribution into long tibble
#'
#' @description The lowerGPreff model functions expect data in a long format,
#'  which is structured to have a value for every unique date and jurisdiction
#'  pair. This function create a tibble of this structure out of a single
#'  distribution that should be replicated in each cell.
#'
#' @param dates infection dates sequence
#' @param jurisdictions jurisdiction names
#' @param value distribution to be replicated in each cell
#'
#' @importFrom dplyr mutate
#' @importFrom tibble tibble
#'
#' @return long tibble with distribution replicated
#' @export
create_lowerGPreff_massfun_timeseries <- function (dates,
                                                jurisdictions,
                                                value) {

  timeseries <- lowerGPreff::create_lowerGPreff_timeseries(
    dates = infection_days,
    jurisdictions = jurisdictions,
    value = list(value))

  class(timeseries) <- c("lowerGPreff_massfun_timeseries",
                         class(value),  # CHECK THIS
                         class(timeseries))
  timeseries

}

#' Title
#'
#' @param dates infection dates sequence
#' @param jurisdictions jurisdiction names
#' @param value constant value to be replicated in each cell
#'
#' @return long tibble with fixed value replicated
#' @export
create_lowerGPreff_fixed_timeseries <- function (dates,
                                                 jurisdictions,
                                                 value) {

  timeseries <- lowerGPreff::create_lowerGPreff_timeseries(
    dates = infection_days,
    jurisdictions = jurisdictions,
    value = value)
  #
  # distribution_timeseries <- dplyr::rename(
  #   notif_full_delay_dist,
  #   distribution = value)

  class(timeseries) <- c("lowerGPreff_fixed_timeseries",
                                      class(timeseries))
  timeseries

}


