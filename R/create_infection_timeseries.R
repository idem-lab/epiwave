#' Create infection timeseries
#'
#' @description Compute the number of infections per day, given a
#'  log-scaled initial number of infections and a time-varying random effect
#'  controlling the trend of infections over time. The random effect is defined
#'  as a Gaussian Process (GP), and can be applied to the infection number
#'  through three different formulations.
#'
#'  This random effect over time is not explicitly specified as an
#'  epidemiologically defined quantity like the growth rate or the reproduction
#'  number. Rather, these quantities are calculated from posterior samples.
#'  The function uses an uninformative prior for the initial number of
#'  infections.
#'
#' @param n_days_infection length of full date sequence of infection timeseries
#' @param n_jurisdictions number of jurisdictions, defaults to 1
#' @param effect_type options include c('infections', 'growth_rate',
#'        'growth_rate_derivative')
#'
#' @importFrom greta lognormal normal
#' @importFrom greta.gp gp mat52
#'
#' @return greta arrays for infection timeseries
#' @export
create_infection_timeseries <- function (n_days_infection,
                                         n_jurisdictions = 1,
                                         effect_type = c('infections',
                                                         'growth_rate',
                                                         'growth_rate_deriv')) {

  # kernel hyperparams
  gp_lengthscale <- greta::lognormal(3, 1)
  gp_variance <- greta::normal(0, 1, truncation = c(0, Inf))
  gp_kernel <- greta.gp::mat52(gp_lengthscale, gp_variance)

  # define gp
  gp <- greta.gp::gp(
    x = seq_len(n_days_infection),
    kernel = gp_kernel,
    n = n_jurisdictions,
    tol = .Machine$double.eps
  )

  # compute infections from gp
  # prior for initial num. of infections on log scale
  inits <- greta::lognormal(0, 1, dim = ncol(gp))

  f <- function (inits, gp, type) {
    switch (type,
            infections = infections(inits, gp),
            growth_rate = growth_rate(inits, gp),
            growth_rate_deriv = growth_rate_deriv(inits, gp))
  }

  infection_timeseries <- f(inits, gp, effect_type)

  greta_arrays <- module(
    gp,
    infection_timeseries,
    gp_lengthscale,
    gp_variance,
    gp_kernel
  )

  return(greta_arrays)
}

#' Effect type "infections" formula for infection timeseries
#'
#' @param inits prior for initial num. of infections on log scale
#' @param z gaussian process
#'
#' @importFrom greta sweep
#'
#' @return infection timeseries
infections <- function (inits, z) {
  exp(greta::sweep(z, 2, inits, FUN = '+'))
}

#' Effect type "growth_rate" formula for infection timeseries
#'
#' @param inits prior for initial num. of infections on log scale
#' @param z gaussian process
#'
#' @importFrom greta apply sweep
#'
#' @return infection timeseries
growth_rate <- function (inits, z) {
  log_rt <- z
  exp(greta::sweep(greta::apply(log_rt, 2, 'cumsum'), 2, inits, FUN = '+'))
}

#' Effect type "growth_rate_deriv" formula for infection timeseries
#'
#' @param inits prior for initial num. of infections on log scale
#' @param z gaussian process
#'
#' @importFrom greta apply sweep
#'
#' @return infection timeseries
growth_rate_deriv <- function (inits, z) {
  log_rt_diff <- z
  log_rt <- greta::apply(log_rt_diff, 2, 'cumsum')
  exp(greta::sweep(greta::apply(log_rt, 2, 'cumsum'), 2, inits, FUN = '+'))
}

