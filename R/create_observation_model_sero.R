#' Create observation model
#'
#' @description This function constructs the observation model and defines
#'   likelihood over the observational data (e.g. cases, hospital
#'   admissions...). It first begets a timeseries of expected observations from
#'   a timeseries of new infections, the infection-to-observation delay
#'   distribution(s), and the proportion of infections to be observed.
#'   Specifically, the model uses the delay distribution(s) to convolve the
#'   infection timeseries into the expected observation timeseries, thus
#'   accounting for the probability of observation over time-since-infection. An
#'   additional multiplication by observation proportion is applied to the
#'   convolved timeseries, to adjust for effects such as case ascertainment.
#'   Optionally, a day-of-week effect may be included in the convolution process
#'   to introduce weekly periodicity. The expected observation timeseries is
#'   treated as the mean of a negative binomial distribution from which the data
#'   is observed, thus allowing us to define likelihood over the data, and link
#'   the data to the unknown infection timeseries.
#'
#' @param sero_decay_curve time-constant decay curve for the sero data used,
#'   expect user to have used GPreff::expand_constant_value to make it the same
#'   time-length as the infection days OR this would be timevarying if
#'   reinfection correction is applied, i.e. the curve at time t takes into
#'   consideration of people infected at time t that would be removed from the
#'   probabiltiy at time delta + t because they are reinfected at time delta + t
#' @param positive_counts counts of seropositive individuals in each survey,
#'   needs by in long format indexed by date and state
#' @param sample_sizes counts of total individuals in each survey, needs by in
#'   long format indexed by date and state
#' @param total_pops total population of the jurisdiction
#' @param infection_timeseries greta array of infection timeseries
#' @param data_id leave this here in case have multiple sero data eg spike vs NP
#'
#' @importFrom greta %*% as_data negative_binomial normal sweep zeros
#'
#' @return greta arrays of observation model
#' @export
create_sero_observation_model <- function (infection_timeseries,
                                           sero_decay_curve,
                                           positive_counts,
                                           sample_sizes,
                                           total_pops,
                                           data_id = NULL
                                           ) {

  # add if statements to check that infection_days is long enough to cover
  # the right period

  positive_counts_mat <- data_to_matrix(positive_counts, 'count')
  decay_mat <- data_to_matrix(sero_decay_curve, 'decay_curve')
  sample_sizes_mat <- data_to_matrix(sample_sizes, 'total')

  infection_days <- as.Date(rownames(decay_mat))

  n_jurisdictions <- ncol(infection_timeseries)
  n_dates <- nrow(infection_timeseries)

  convolution_matrices <- lapply(
    1:n_jurisdictions,
    function(x) {
      get_convolution_matrix(decay_mat[, x],
                             n_dates)
    })

  # compute expected seropositives of the same length
  expected_positive_list <- lapply(
    1:n_jurisdictions,
    function(x) {
      convolution_matrices[[x]] %*% infection_timeseries[, x]
    })

  expected_positive <- do.call(
    cbind,
    expected_positive_list
  )

  # divide by total pop of jurisdictions to get pop positivity rates
  expected_seropositivity_rate <- sweep(expected_positive,
                                        2,
                                        total_pops,
                                        FUN = "/")

  # get time idx for when sampling took place
  sampling_idx <- infection_days %in% rownames(positive_counts_mat)
  # subset to expected positives for the sampled days
  expected_seropositivity_rate_sampled <- expected_seropositivity_rate[sampling_idx, ]

  valid_mat <- positive_counts_mat
  valid_mat[is.na(positive_counts_mat)] <- FALSE
  valid_mat[!is.na(positive_counts_mat)] <- TRUE
  valid_idx <- as.logical(as.numeric(valid_mat))

  positive_counts_mat_array <- greta::as_data(
    as.numeric(positive_counts_mat)[valid_idx])

  greta::distribution(positive_counts_mat_array) <- greta::binomial(
    size = sample_sizes_mat[valid_idx],
    prob = expected_seropositivity_rate_sampled[valid_idx]
  )

  greta_arrays <- list(
    expected_seropositivity_rate, # suspect we want to output the full length version?
    convolution_matrices
  )

  names(greta_arrays) <- c(
    paste0(data_id, '_expected_seropositivity_rate'),
    paste0(data_id, '_convolution_matrices')
  )

  return(greta_arrays)
}
