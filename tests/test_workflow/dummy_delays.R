dummy_delay1 <- data.frame(delays = c(1:4),
                          mass = c(0.25, 0.5, 0.25, 0))
class(dummy_delay1) <- c("epiwave_distribution_massfun", class(dummy_delay1))

dummy_delay2 <- data.frame(delays = c(1:3),
                           mass = c(0.25, 0.5, 0.25))
class(dummy_delay2) <- c("epiwave_distribution_massfun", class(dummy_delay2))


# amend plot so that days without inference aren't shown -- fill with NAs.
