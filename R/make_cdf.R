#' Make cumulative density function
#'
#' @description Define cumulative density function. Used especially for
#'  calculation of incubation period. Parameters for Omicron and Delta
#'  SARS-CoV-2 strains estimated from
#'  https://www.eurosurveillance.org/content/10.2807/1560-7917.ES.2022.27.6.2200042.
#'  This is not a time-varying biological quantity - any changes in the
#'  dominant variant should be modelled separately.
#'
#' @param option SARS-CoV-2 covariant or NULL model
#'
#' @importFrom stats approxfun pweibull
#'
#' @return cumulative density function for incubation period
#' @export
make_cdf <- function (option = c('Omicron',
                                 'Delta',
                                 'None')) {

  days <- 0:28

  f <- function (option) {
    switch (option,
            Delta = stats::pweibull(days, shape = 1.83, scale = 4.93),
            Omicron = stats::pweibull(days, shape = 1.5, scale = 3.6),
            None = stats::pweibull(days, shape = 1, scale = 0.1))
  }

  cum_density <- f(option)
  cdf <- stats::approxfun(days, cum_density)

  return(cdf)
}
