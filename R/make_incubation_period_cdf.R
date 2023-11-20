#' Make incubation period cumulative density function
#'
#' @description Note that in using incubation period for modelling we should
#' not treat it as time-varying (due to changes in the dominant variant).
#' Maybe through vignette or some other documentation we should make it clear
#' that this is not a time-varying biological quantity and the change in
#' dominant variant should be modelled separately)
#'
#' @param strain
#'
#' @importFrom stats approxfun pweibull
#'
#' @return cumulative density function for incubation period
#' @export
make_incubation_period_cdf <- function (strain = c('Omicron',
                                                   'WT',
                                                   'Alpha',
                                                   'Beta/Gamma',
                                                   'Delta')) {

  strain <- match.arg(strain)

  if (strain == 'Omicron') {

    # parameters estimated from
    # https://www.eurosurveillance.org/content/10.2807/1560-7917.ES.2022.27.6.2200042
    days <- 0:28
    cum_density <- stats::pweibull(days, shape = 1.5, scale = 3.6)

  }

  cdf <- stats::approxfun(days, cum_density)

  return(cdf)
}
