# this series of related functions help construct discrete probability mass
# functions for time intervals, which are used for the convolution of
# timeserieses. For a range of real-valued integer time t, a known timeseries
# A_t for event A and a time interval mass function G(t',t), we want to get a
# convolved timeseries for event B, B_t, such that B_t = Sum_(t'=0)^t
# G(t',t)*A_t', that is, the value of B_t is the sum of the values of A at all
# time points prior to and including t, multiplied by the probability of the A
# to B event delay having the length of t-t'. This probability is given by the
# convolution mass function G(t',t), and can vary with t. This assumes that the
# event A occurs no later than event B, e.g convolving from count timeseries of
# infection to count timeseries of incidence notification,  or from count
# timeseries of infection to pop-level seropositivity timeseries. If event B can
# predate event A (e.g. notification before symptom), then the sum needs to be
# computed over Sum_(t'=0)^(t+x) with x being the plausible minimum of the
# negative delay value. Note this requires the convolution mass function to have
# some mass over negative delay-time values. Different versions of these
# functions enable different types of input information. Principally, they
# differ by if the input information enables the construction of a parametric
# distribution (e.g. epistudies proposing parametric distributions for the
# incubation period or time delay from infection to hospitalisation for
# hospitalised cases), or if empirical time mass is observable, but no
# parametric shapes can be easily proposed. In the latter case, we build a
# discrete empirical mass function, either directly from the observed time
# interval data, or convert from a pre-calculated empirical cdf. The cdf input
# is considered here because it allows one to share this information without
# sharing other potentially sensitive aspects of the data from which the time
# intervals are observed from. Note that for some situation such as convolving
# from infection count timeseries to notification count timeseries, we can
# assume that every infection must eventually lead to a notification (as we are
# dealing with ascertainment elsewhere), therefore G(t',tau) is a strict
# probability mass function and sums to 1 over a reasonably large range of tau,
# and Sum_t B_t is roughly equal to Sum_t A_t save for the right truncation
# caveat. In other cases, such as convolving from infections to pop-level
# seropositivity, the infected can test seropositive multiple times, so
# G(t',tau) is not a probability mass function, but rather describes the decay
# in seropositivity rate. That is, the seropositivity curve describes the
# probability of testing positive at time tau-t' post infection. This has
# implications on the identifiability of B_t.


# need to think about parametric form, should need a separate function for each distribution?
conv_mass_weibull <- function(min_delay = 0,
                              max_delay = 56, #need to think about this carefully, this could be a strong assumption skewing and truncating the shape of the delay
                              output_is_pmf = TRUE,
                              shape = NULL,
                              scale = NULL) {
  x <- min_delay:max_delay
  y <- stats::dweibull(x,
                       shape,
                       scale)
  # if the dist is supposed to be a pmf, make sure it sums to 1
  if (output_is_pmf) {
    y <- y/sum(y)
  }

  stats::approxfun(x, y)
}

# this is a simple conversion from cdf
conv_mass_from_cdf <- function(min_delay = 0,
                               max_delay = 56,
                               output_is_pmf = TRUE,
                               input_cdf){

  x <- min_delay:max_delay
  # calculate prob at each time step
  y <- input_cdf(x + 1) - input_cdf(x)
  # if the dist is supposed to be a pmf, make sure it sums to 1
  if (is_pmf) {
    y <- y/sum(y)
  }
  stats::approxfun(x, y)
}

# calculate from empirical delay data
conv_mass_empirical <- function(min_delay = -3, #example for when anticipating possible negative delay
                                max_delay = 56,
                                input_data){

  x <- min_delay:max_delay

  if (!is.numeric(input_data)) stop("input data should be time intervals (e.g. delays) but don't look like numbers!")

  input_data <- na.omit(input_data)
  # calculate cdf from the data
  emp_cdf <- stats::ecdf(input_data)
  # calculate prob at each time step
  y <- emp_cdf(x + 1) - emp_cdf(x)
  # if the dist is supposed to be a pmf, make sure it sums to 1
  if (output_is_pmf) {
    y <- y/sum(y)
  }
  stats::approxfun(x, y)
}
