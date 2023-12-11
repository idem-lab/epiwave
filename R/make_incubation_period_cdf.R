#' Make incubation period cumulative density function
#'
#' @description Define cumulative density function of incubation period.
#'  Parameters for different COVID-19 strains estimated from
#'  https://www.eurosurveillance.org/content/10.2807/1560-7917.ES.2022.27.6.2200042.
#'  This is not a time-varying biological quantity - any changes in the
#'  dominant variant should be modelled separately.
#'
#' @param strain COVID-19 strain
#'
#' @importFrom stats approxfun pweibull
#'
#' @return cumulative density function for incubation period
#' @export
make_incubation_period_cdf <- function (strain = c('Omicron',
                                                   # 'WT',
                                                   # 'Alpha',
                                                   # 'Beta/Gamma',
                                                   'Delta')) {

  strain <- match.arg(strain)

  if (strain == 'Omicron') {

    days <- 0:28
    cum_density <- stats::pweibull(days, shape = 1.5, scale = 3.6)

  }

  if (strain == 'Delta') {

  }

  cdf <- stats::approxfun(days, cum_density)

  return(cdf)
}
