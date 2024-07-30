fit_waves <- function (observations, # list of data type lists
                      target_infection_dates = NULL) {

  # prep the model objects

  # delays <- observations$delays

  # add check that ncol(infection_timeseries and below yield same. number of juris)

  n_jurisdictions <- 8#length(unique(delays$jurisdiction))
  n_days_infection <- length(target_infection_dates)



  # sanitise dates
  # target_infection_dates <- check_target_infection_dates(target_infection_dates)
  # observations <- sanitise_dates(observations, target_infection_dates)


  # set up the greta model
  ## infection timeseries model
  infection_model_objects <- epiwave::create_infection_timeseries(
    n_days_infection,
    n_jurisdictions,
    effect_type = 'growth_rate')

  # observtion model objects in observations
  test <- lapply(names(observations),
                 create_observation_model,
                 observations, infection_model_objects)

  # observation_model_objects <- new_create_observation_model(
  #   observations,
  #   infection_model_objects,
  #   data_id = 'cases')


  ## greta model fit
  m <- greta::model(infection_model_objects$infection_timeseries)

  fit <- epiwave.pipelines::fit_model(
    model = m,
    n_chains = 2,
    max_convergence_tries = 1,
    warmup = 100,
    n_samples = 100,
    n_extra_samples = 100)

  # return the stuff in an object of type waveylistythingy
  wavylistythingy <- list(
    infection_model_objects = infection_model_objects,
    fit = fit,
    infection_days = infection_days,
    jurisdictions = jurisdictions)

  return(wavylistythingy)

}
