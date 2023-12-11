#' Implement day of week effect in observation model
#'
#' @param infection_days vector of dates that cover infection timeseries
#' @param dow_model greta model for day of week effect
#'
#' @return matrix of dimension date x jurisdiction with day of week correction
#' @export
implement_day_of_week <- function (infection_days,
                                   dow_model) {

  # get day of week index
  # favour base approach over lubridate for manual declaration of factor level
  # so that we know which index is which day of the week
  dweek <- weekdays(infection_days)
  dweek <- as.integer(
    factor(
      dweek,
      levels = c("Monday",
                 "Tuesday",
                 "Wednesday",
                 "Thursday",
                 "Friday",
                 "Saturday",
                 "Sunday")
    )
  )

  # normalise multiplier to average to 1
  dow_weights <- dow_model$dow_dist * 7

  # match weight to date by jurisdiction matrix
  dow_correction <- t(dow_weights[, dweek])

  return(dow_correction)
}
