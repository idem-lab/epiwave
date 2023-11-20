#' Make generation interval density
#'
#' @param gi_distribution_data
#'
#' @importFrom stats plnorm
#'
#' @return
#' @export
make_generation_interval_density <- function (gi_distribution_data) {

    # generation interval distribution; use SI distribution from Nishiura et al.
    meanlog <- mean(gi_distribution_data$param1)
    sdlog <- mean(gi_distribution_data$param2)

    # this has to be kept inside the parent function so the meanlog and
    # sdlog pass correctly
    gi_density <- function (days) {
        stats::plnorm(days, meanlog, sdlog)
    }

    gi_out <- construct_delays(gi_density,
                               output = "probability",
                               stepfun_output = TRUE)

    gi_out
}

