#' Create infection timeseries
#'
#' @description Compute the number of infections per day either with an
#'  uninformative flat prior (effect_type = 'flat_prior'), or given a
#'  log-scaled initial number of infections and a time-varying random effect
#'  controlling the trend of infections over time. The random effect is defined
#'  as a Gaussian Process (GP), and can be applied to the infection number
#'  through three different formulations: effect_type = 'gp_infections',
#'  'gp_growth_rate', or 'gp_growth_rate_derivative'.
#'
#'  gp_infections: the log of infections tend towards the GP mean value in the
#'  now-cast period.
#'  gp_growth_rate: the growth rate tends to 1 and the infections stay at
#'  around the same level in now-casts and forecasts.
#'  gp_growth_rate_derivative: in now-cast/forecast the infection trajectory
#'  would follow the most recent growth rate trend.
#'
#'  This random effect over time is not explicitly specified as an
#'  epidemiologically defined quantity like the growth rate or the reproduction
#'  number. Rather, these quantities are calculated from posterior samples.
#'  The function uses an uninformative prior for the initial number of
#'  infections.
#'
#' @param n_days_infection length of full date sequence of infection timeseries
#' @param n_jurisdictions number of jurisdictions, defaults to 1
#' @param effect_type options include 'flat_prior', 'gp_infections', 'gp_growth_rate',
#'        'gp_growth_rate_derivative'. See description for more info.
#'
#' @importFrom greta lognormal normal
#' @importFrom greta.gp gp mat52
#'
#' @return greta arrays for infection timeseries
#' @export
create_infection_timeseries <- function (n_days_infection,
                                         n_jurisdictions = 1,
                                         effect_type = c('flat_prior',
                                                         'gp_infections',
                                                         'gp_growth_rate',
                                                         'gp_growth_rate_deriv')) {

  # improper flat prior for infection incidence
  if (effect_type == 'flat_prior') {
    infection_timeseries <- greta::variable(lower = 0,
                                            dim = c(n_days_infection,
                                                    n_jurisdictions))
    greta_arrays <- module(
      infection_timeseries
    )
    return(greta_arrays)
  }

  # kernel hyperparams
  gp_lengthscale <- greta::lognormal(0, 3) #inverse_gamma(187/9,1157/18)
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
            gp_infections = gp_infections(inits, gp),
            gp_growth_rate = gp_growth_rate(inits, gp),
            gp_growth_rate_deriv = gp_growth_rate_deriv(inits, gp))
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
gp_infections <- function (inits, z) {
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
gp_growth_rate <- function (inits, z) {
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
gp_growth_rate_deriv <- function (inits, z) {
  log_rt_diff <- z
  log_rt <- greta::apply(log_rt_diff, 2, 'cumsum')
  exp(greta::sweep(greta::apply(log_rt, 2, 'cumsum'), 2, inits, FUN = '+'))
}

