#' Title
#'
#' @param delay1
#' @param delay2
#'
#' @return
#' @export
#'
#' @examples
add_distributions <- function (delay1, delay2) {

  added_dists <- add(delay1, delay2)

  timeseries_added_dists <- create_epiwave_massfun_timeseries(
    dates = infection_days,
    jurisdictions = jurisdictions,
    value = added_dists)
  timeseries_added_dists

}
