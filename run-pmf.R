source('R/pmf_dev.R')
library(dplyr)

dist <- distributional::dist_weibull(shape = 2, scale = 3)
delay_probs_from_dist <- params_as_delay_dist(dist)



data <- data.frame(sym_date = rep(seq(from = as.Date('2023-01-2'),
                                      to = as.Date('2023-03-4'),
                                      by = 'day'),
                                  each = 20))
data$notif_date <- data$sym_date + rpois(nrow(data), 4)
delay_probs_from_data <- data_as_delay_dist(data)


combined <- combine_massfuns(delay_probs_from_dist, delay_probs_from_data)
