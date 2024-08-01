fit_waves <- function (observations, # list of data type lists
                      target_infection_dates = NULL,
                      n_chains = 2,
                      max_convergence_tries = 1,
                      warmup = 100,
                      n_samples = 100,
                      n_extra_samples = 100) {

  # prep the model objects

  # delays <- observations$delays

  # add check that ncol(infection_timeseries and below yield same. number of juris)


  n_jurisdictions <- length(unique(observations[[1]]$timeseries_data$jurisdiction))
  n_days_infection <- length(target_infection_dates)



  # sanitise dates
  # target_infection_dates <- check_target_infection_dates(target_infection_dates)
  # observations <- sanitise_dates(observations, target_infection_dates)


  # set up the greta model
  ## infection timeseries model
  infection_model <- epiwave::create_infection_timeseries(
    n_days_infection,
    n_jurisdictions,
    effect_type = 'growth_rate')

  # observtion model objects in observations
  observation_models <- lapply(names(observations),
                               create_observation_model,
                               observations, infection_model)

  # case_model_objects <- create_observation_model(
  #   names(observations)[1],
  #   observations,
  #   infection_model)
  # hosp_model_objects <- create_observation_model(
  #   names(observations)[2],
  #   observations,
  #   infection_model)




  ## greta model fit
  m <- greta::model(infection_model$infection_timeseries)

  fit <- epiwave.pipelines::fit_model(
    model = m,
    n_chains = n_chains,
    max_convergence_tries = max_convergence_tries,
    warmup = warmup,
    n_samples = n_samples,
    n_extra_samples = n_extra_samples)

  # return the stuff in an object of type waveylistythingy
  wavylistythingy <- list(
    infection_model_objects = infection_model_objects,
    fit = fit,
    infection_days = infection_days,
    jurisdictions = jurisdictions)

  return(wavylistythingy)

}
