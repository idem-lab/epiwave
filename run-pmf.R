# source('R/pmf_dev.R')
library(dplyr)
library(lowerGPreff)

dist <- distributional::dist_weibull(shape = 2, scale = 3)
delay_probs_from_dist <- parametric_dist_to_distribution(dist)



data <- data.frame(sym_date = rep(seq(from = as.Date('2023-01-2'),
                                      to = as.Date('2023-03-4'),
                                      by = 'day'),
                                  each = 20))
data$notif_date <- data$sym_date + rpois(nrow(data), 4)
delay_probs_from_data <- data_to_distribution(data)


combined <- combine(delay_probs_from_dist,
                    delay_probs_from_data)
