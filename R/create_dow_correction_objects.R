#' Create day of week correction arrays
#'
#' @param infection_dates full date sequence of infection timeseries
#' @param n_jurisdictions number of jurisdictions, defaults to 1
#' @param data_id optional name label identifying data type for greta arrays
#'
#' @return named greta arrays for day of week effects
#' @export
create_dow_correction_arrays <- function(
        infection_dates,
        n_jurisdictions = 1,
        data_id = '') {

    # get day of week index
    # favour base approach over lubridate for manual declaration of factor level
    # so that we know which index is which day of the week
    dweek <- weekdays(infection_dates)
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

    # prior for dweek correction
    dow_alpha <- greta::normal(1, 1,
                               truncation = c(0, Inf),
                               dim = c(1, 7))

    dow_dist <- greta::dirichlet(dow_alpha,
                                 n_realisations = n_jurisdictions)

    # normalise multiplier to average to 1
    dow_weights <- dow_dist * 7

    # match weight to date by state matrix
    dow_correction <- t(dow_weights[, dweek])

    greta_arrays <- list(
        dow_weights,
        dow_correction)

    names(greta_arrays) <- c(
        paste0(data_id, '_dow_weights'),
        paste0(data_id, '_dow_correction'))

    return(greta_arrays)
}
