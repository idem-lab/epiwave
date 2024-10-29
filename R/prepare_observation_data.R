#' Prepare observation data
#'
#' @description Prepare list of data objects needed for observation model.
#'
#' @param timeseries_data long form data of counts of notifications
#' @param delay_from_infection distribution of delays, covering length of
#'  infection timeseries, including incubation period, if applicable
#' @param proportion_infections single fixed value of expected proportions of
#'   infections observed in count data
#' @param type count or prevalence
#' @param dow_model optional module of greta arrays defining day-of-week model
#'
#' @importFrom greta %*% as_data negative_binomial normal sweep zeros
#'
#' @return greta arrays of observation model
#' @export
prepare_observation_data <- function (timeseries_data,
                                      delay_from_infection,
                                      proportion_infections,
                                      type = c('count', 'prevalence'),
                                      dow_model = NULL) {

  # see how we are defining the likelihood
  type <- match.arg(type)


  # check the data is a valid timeseries object (or can be coerced to one) and
  # that it contains data consistent with the observation model type
  # timeseries_data <- check_valid_timeseries_data(timeseries_data,
  #                                                type = type)

  # add if statements to check that infection_days is long enough to cover
  # the right period

  # these should be two different dimensions yes? because prop will be
  # dim of infection dates and cases will be cases.
  # if these should both be infection dates we can change things.
  # because this doesn't work for a single value for prop.
  # even though we can pass a single value through to the convolution
  # we can do the cheat below to get back infection dates for dow.
  case_mat <- as_matrix(timeseries_data)
  prop_mat <- as_matrix(proportion_infections)


  # is this basically a cheat way to get back the infection dates?
  # if so it's sloppy, if not, are these supposed to match case dates?
  obs_infection_days <- as.Date(rownames(prop_mat))

  ## add a check for correct dow arrays
  if (!is.null(dow_model)) {
    dow_correction <- implement_day_of_week(obs_infection_days, dow_model)
    prop_mat <- prop_mat * dow_correction
  }

  out <- list(timeseries_data = timeseries_data,
              delays = delay_from_infection,
              case_mat = case_mat,
              prop_mat = prop_mat)
  out

}
