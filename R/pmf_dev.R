match_mass_lookup <- function (day_diff, lookup, mass_val_col, delay_col) {

  lookup <- as.data.frame(lookup)
  day_diff[] <- lookup[,mass_val_col][
    match(unlist(day_diff),
          lookup[,delay_col])]
  day_diff[is.na(day_diff)] <- 0

  class(day_diff) <- 'numeric'
  day_diff

}

make_massfun <- function (min_delay,
                          max_delay,
                          cdf_fun,
                          normalise = c(TRUE, FALSE)) {

  delay_massfun <- data.frame(
    delays = seq(min_delay, max_delay)
  ) %>%
    dplyr::mutate(
      upper = cdf_fun(delays + 1),
      lower = cdf_fun(delays),
      mass = upper - lower
    )

  if (normalise) {
    delay_massfun <- delay_massfun %>%
      dplyr::mutate(
        correction = sum(mass),
        mass = mass / correction
      )
  }

  delay_massfun <- delay_massfun %>%
    dplyr::select(
      delays, mass
    )

  # add if to separate the two classes that are normlaise or not.
  class(delay_massfun) <- c("delay_probs", class(delay_massfun))

  delay_massfun

}

# potential new name = data_to_distribution
data_as_delay_dist <- function (data,
                                min_delay = NULL,
                                max_delay = NULL,
                                normalise = TRUE) {

  day_diff <- data$notif_date - data$sym_date
  cdf_fun <- stats::ecdf(day_diff)

  if (is.null(min_delay)) {
    min_delay <- floor(quantile(cdf_fun, 0.00))
  }
  if (is.null(max_delay)) {
    max_delay <- ceiling(quantile(cdf_fun, 0.99))
  }

  out <- make_massfun(min_delay, max_delay, cdf_fun, normalise = normalise)
  out

}

# potential new name = parametric_dist_to_distribution
params_as_delay_dist <- function (dist,
                                  min_delay = NULL,
                                  max_delay = NULL,
                                  normalise = TRUE) {

  # this quantile is the S3 method for quantile from distributional pkg
  if (is.null(min_delay)) {
    min_delay <- floor(quantile(dist, 0.00))
  }
  if (is.null(max_delay)) {
    max_delay <- ceiling(quantile(dist, 0.99))
  }

  # when distributional::cdf is applied to a sequence (of x) it returns a list
  cdf_fun <- function(x) distributional::cdf(dist, x)[[1]]

  out <- make_massfun(min_delay, max_delay, cdf_fun, normalise = normalise)
  out

}

combine_massfuns <- function (delay_massfun1, delay_massfun2) {

  names(delay_massfun1) <- c('delay1', 'massfun1')
  names(delay_massfun2) <- c('delay2', 'massfun2')

  p <- tidyr::expand_grid(
    delay1 = delay_massfun1$delay1, # ie. incubation period
    delay2 = delay_massfun2$delay2) %>% # ie. sym to notif
    dplyr::left_join(delay_massfun1) %>%
    dplyr::left_join(delay_massfun2) %>%
    dplyr::mutate(total_delay = delay1 + delay2) %>%
    dplyr::group_by(total_delay) %>%
    dplyr::summarise(massfun_combined =
                       sum(massfun1 * massfun2))

  p

}


plot.delay_probs <- function(x, y, ...) {
  barplot(x$mass, width = 1, names.arg = x$delays,
          xlab = "delay (days)",
          ylab = "probability")
}
