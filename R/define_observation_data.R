define_observation_data <- function (timeseries_data,
                                     delay_from_infection,
                                     proportion_infections,
                                     dow_model = NULL) {



  out <- list(timeseries_data = timeseries_data,
              delay_from_infection = delay_from_infection,
              proportion_infections = proportion_infections,
              dow_model = dow_model)
  out

}
