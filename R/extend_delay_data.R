#' Extend delay data
#'
#' @param delay_dat long form delay data that covers the date sequence of the
#'  entire infection timeseries
#' @param incubation_period_distribution distribution for the incubation
#'  period
#'
#' @return delay distribution that combines both infection to notification
#'  delay and incubation period
#'
#' @export
extend_delay_data <- function (delay_dat,
                               incubation_period_distribution) {

  delay_matrix <- data_to_matrix(delay_dat, 'delay_fxn')

  delay_and_incubation <- apply(
    delay_matrix,
    c(1,2), construct_delays,
    ecdf2 = incubation_period_distribution,
    output = "probability",
    stepfun_output = TRUE)

  # make so this function spits back out long. so change the construct_delays
  # to just add the two delays.

  delay_and_incubation

  #### if the delay_dat input is the right dates below is not needed.

  # earliest_date <- as.Date(min(rownames(delay_and_incubation)))
  # latest_date <- as.Date(max(rownames(delay_and_incubation)))
  #
  # extra_left <- as.character(earliest_date - min(infection_days))
  # extra_right <- as.character(max(infection_days) - latest_date)
  #
  # delay_matrix_ext_left <- rbind(
  #   do.call("rbind",
  #           replicate(extra_left,
  #                     delay_and_incubation[1,],
  #                     simplify = FALSE)),
  #   delay_and_incubation)
  #
  # delay_matrix_ext_both <- rbind(
  #   delay_matrix_ext_left,
  #   do.call("rbind",
  #           replicate(extra_right,
  #                     delay_and_incubation[nrow(delay_and_incubation),],
  #                     simplify = FALSE)))
  #
  # rownames(delay_matrix_ext_both) <- as.character(infection_days)
  #
  # delay_matrix_ext_both
}

