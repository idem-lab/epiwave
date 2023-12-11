#' @description  Make generation interval density following data from the
#'   literature. This function requires inputs of mean and sd estimates of the
#'   generation interval on the log scale, and outputs a GI distribution in the
#'   bespoke pmf format used by this package. By default this function expects
#'   input relevant to COVID-19, specifically, samples provided by Nishiura et
#'   al 2020 (https://doi.org/10.1016/j.ijid.2020.02.060). This estimation is
#'   for serial interval, not generation interval, since the latter was
#'   difficult to observe at the time of the Nishiura publication. However,
#'   comparison against more recent household studies shows the Nishiura
#'   estimates still accurrately reflect our current understanding of likely
#'   COVID GI.
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

