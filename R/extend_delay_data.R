#' Extend delay data
#'
#' @param delay_dat long form delay data that covers the date sequence of the
#'  entire infection timeseries
#' @param incubation_period_distribution optional distribution for the
#'  incubation period
#'
#' @return delay distribution that combines both infection to notification
#'  delay and incubation period
#'
#' @export
extend_delay_data <- function (delay_dat,
                               incubation_period_distribution = NULL) {

  delay_matrix <- data_to_matrix(delay_dat, 'delay_fxn')

  # our short-term solution to allow hospitalisation data to be reformatted
  # with this function, even though doesn't need to be merged with incubation
  if (is.null(incubation_period_distribution)) {
    incubation_period_distribution <- make_cdf(option = "None")
  }

  delay_and_incubation <- apply(
    delay_matrix,
    c(1,2), construct_delays,
    ecdf2 = incubation_period_distribution,
    output = "probability",
    stepfun_output = TRUE)

  delay_and_incubation
}

