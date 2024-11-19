#' Fit waves to observations
#'
#' @description
#'
#' Compute the number of infections per day either with an
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
#' @param observations prepared datasets
#' @param infection_model_type options include 'flat_prior', 'infections', 'growth_rate',
#'        'growth_rate_derivative'. See description for more info.
#' @param target_infection_dates infection dates that cover more than the data dates
#' @param n_chains number of chains to run MCMC
#' @param max_convergence_tries number of times to repeat run to get convergence
#' @param warmup how much samples to use in warmup
#' @param n_samples number of samples after warmup
#' @param n_extra_samples number of extra samples if original run didn't converge
#'
#' @return list of infection_model, fit, infection_days, and jurisdictions
#' @export
#'
fit_waves <- function (observations,
                       infection_model_type = c('flat_prior',
                                                'gp_infections',
                                                'gp_growth_rate',
                                                'gp_growth_rate_deriv'),
                       n_chains = 4,
                       max_convergence_tries = 3,
                       warmup = 1000,
                       n_samples = 2000,
                       n_extra_samples = 1000) {


  # prep the model objects
  # add check that ncol(infection_timeseries and below yield same. number of juris)

  target_infection_dates <- observations$target_infection_dates
  n_days_infection <- length(target_infection_dates)

  jurisdictions <- observations$target_jurisdictions
  n_jurisdictions <- length(jurisdictions)

  # set up the greta model
  # infection timeseries model
  incidence_greta_arrays <- create_infection_timeseries(
    n_days_infection,
    n_jurisdictions,
    effect_type = infection_model_type)
  incidence <- incidence_greta_arrays$infection_timeseries

  # observation model objects in observations
  observation_model_data <- observations$observation_model_data
  observation_models <- lapply(names(observation_model_data),
                               create_observation_model,
                               observation_model_data,
                               target_infection_dates,
                               incidence)
  names(observation_models) <- names(observation_model_data)

  # greta model fit
  m <- greta::model(incidence)

  if (infection_model_type == 'flat_prior') {

    inits_df <- observations$inits_df
    first_init <- inits_df[1]
    extra_beginning <- round((n_days_infection - length(inits_df)) / 2)
    last_init <- inits_df[length(inits_df)]
    extra_end <- n_days_infection - (extra_beginning + length(inits_df))
    inits_full <- c(rep(first_init, times = extra_beginning),
                    inits_df,
                    rep(last_init, times = extra_end))

    inits <- greta::initials(
      incidence = inits_full)

  } else { inits <- greta::initials() }

  fit <- greta::mcmc(
    m,
    sampler = greta::hmc(Lmin = 25, Lmax = 30),
    chains = n_chains,
    warmup = warmup,
    n_samples = n_samples,
    initial_values = inits,
    one_by_one = TRUE
  )

  # return the outputs
  fit_output <- list(
    infection_model = incidence,
    fit = fit,
    infection_days = target_infection_dates,
    jurisdictions = jurisdictions)

  return(fit_output)

}
