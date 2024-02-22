#' Expand constant value into long tibble
#'
#' @description The lowerGPreff model functions expect data in a long format,
#'  which is structed to have a value for every unique date and jurisdiction
#'  pair. This function create a tibble of this structure out of a single
#'  value that should be replicated in each cell.
#'
#' @param dates infection dates sequence
#' @param jurisdictions jurisdiction names
#' @param constant_val value to be replicated in each cell
#' @param col_name column name to contain the `constant value`
#'
#' @importFrom dplyr mutate
#' @importFrom tibble tibble
#'
#' @return long tibble with constant value replicated
#' @export
#Assign superclass structure for timeseries data (i.e. distributions)
create_lowerGPreff_timeseries <- function (dates,
                                   jurisdictions,
                                   constant_val) {

  long_unique <- expand.grid(date = dates,
                             jurisdiction = jurisdictions)
  long_combined <- long_unique |>
    dplyr::mutate(value = constant_val) |> #consider another name for constant_val
    #dplyr::mutate(!!col_name := constant_val) |>
    tibble::tibble()

  class(long_combined) <- c("lowerGPreff_timeseries", class(long_combined)) #adds superclass argument

  long_combined

  #also need another function that creates a proper time varying time series (i.e. not outputting constant value)
  #:= name injection, assigns variable name to column name
}


create_lowerGPreff_dist_timeseries <- function (

){
    notif_full_delay_dist <- lowerGPreff::create_lowerGPreff_timeseries(
    dates = infection_days,
    jurisdictions = jurisdictions,
    constant_val = list(notif_delay_dist))

  distribution_timeseries <- dplyr::rename(notif_full_delay_dist, distribution = value)
  class(distribution_timeseries) <- c("lowerGPreff_distribution_timeseries", class(distribution_timeseries)) #adds next level class argument
  distribution_timeseries
}





