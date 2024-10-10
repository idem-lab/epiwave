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
                       target_infection_dates = NULL,
                       n_chains = 4,
                       max_convergence_tries = 3,
                       warmup = 1000,
                       n_samples = 2000,
                       n_extra_samples = 1000) {

  # prep the model objects
  # add check that ncol(infection_timeseries and below yield same. number of juris)

  jurisdictions <- unique(observations[[1]]$timeseries_data$jurisdiction)
  n_jurisdictions <- length(jurisdictions)
  n_days_infection <- length(target_infection_dates)

  # sanitise dates
  # target_infection_dates <- check_target_infection_dates(target_infection_dates)
  # observations <- sanitise_dates(observations, target_infection_dates)

  # set up the greta model
  # infection timeseries model
  infection_model <- epiwave::create_infection_timeseries(
    n_days_infection,
    n_jurisdictions,
    effect_type = infection_model_type)

  # observation model objects in observations
  observation_models <- lapply(names(observations),
                               create_observation_model,
                               observations,
                               target_infection_dates,
                               infection_model)

  # greta model fit
  m <- greta::model(infection_model$infection_timeseries)

  if (infection_model_type == 'flat_prior') one_by_one <- FALSE
  else one_by_one <- TRUE

  fit <- epiwave.pipelines::fit_model(
    model = m,
    n_chains = n_chains,
    max_convergence_tries = max_convergence_tries,
    warmup = warmup,
    n_samples = n_samples,
    n_extra_samples = n_extra_samples,
    one_by_one = one_by_one)

  # return the outputs
  fit_output <- list(
    infection_model = infection_model,
    fit = fit,
    infection_days = target_infection_dates,
    jurisdictions = jurisdictions)

  return(fit_output)

}
