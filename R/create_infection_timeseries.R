#' Create infection timeseries
#'
#' @description WIP function to compute the # of infections per day, given a
#' log-scaled initial # of infections and a time-varying random effect
#' controlling the trend of infections over time. The random effect is defined
#' as a Gaussian Process (GP), and can be applied to the infection # through
#' three different formulations. (see vignette for details here or something).
#' Note that this random effect over time is not explicitly specified as an
#' epidemiologically defined quantity like the growth rate or the reproduction
#' number. Rather, these quantities are to be calculated from posterior samples.
#' The function uses an uninformative prior for the initial # of infections ---
#' currently there are no plans to make this an input param.
#'
#' @param n_days_infection length of full date sequence of infection timeseries
#' @param n_jurisdictions number of jurisdictions, defaults to 1
#' @param effect_type options include c('infections', 'growth_rate',
#'        'growth_rate_derivative')
#'
#' @return greta arrays for infection timeseries
#' @export
create_infection_timeseries <- function (n_days_infection,
                                         n_jurisdictions = 1,
                                         effect_type = c('infections',
                                                         'growth_rate',
                                                         'growth_rate_derivative')) {

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
  inits <- lognormal(0, 1, dim = ncol(gp))

  f <- function (inits, gp, type) {
    switch (type,
            infections = infections(inits, gp),
            growth_rate = growth_rate(inits, gp),
            growth_rate_derivative = growth_rate_derivative(inits, gp))
  }

  infection_timeseries <- f(inits, gp, effect_type)

  greta_arrays <- module(
    gp,
    infection_timeseries, # greta_array but with date/state attributes
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
#' @return infection timeseries
infections <- function (inits, z) {
  exp(sweep(z, 2, inits, FUN = '+'))
}

#' Effect type "growth_rate" formula for infection timeseries
#'
#' @param inits prior for initial num. of infections on log scale
#' @param z gaussian process
#'
#' @return infection timeseries
growth_rate <- function (inits, z) {
  log_rt <- z
  exp(sweep(greta::apply(log_rt, 2, 'cumsum'), 2, inits, FUN = '+'))
}

#' Effect type "growth_rate_derivative" formula for infection timeseries
#'
#' @param inits prior for initial num. of infections on log scale
#' @param z gaussian process
#'
#' @return infection timeseries
growth_rate_derivative <- function (inits, z) {
  log_rt_diff <- z
  log_rt <- greta::apply(log_rt_diff, 2, 'cumsum')
  exp(sweep(greta::apply(log_rt, 2, 'cumsum'), 2, inits, FUN = '+'))
}

