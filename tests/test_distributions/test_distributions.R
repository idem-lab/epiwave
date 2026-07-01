## incubation period
library(epiwave.params)
not_synced_folder <- '../data'

incub_dist <- distributional::dist_weibull(shape = 1.83, scale = 4.93)
incubation <- as_discrete_pmf(incub_dist)
saveRDS(incubation, 'tests/test_distributions/incubation.rds')

## generation interval distribution
gi_distribution_data <- readr::read_csv(
  file = paste0(not_synced_folder,'/nishiura_samples.csv'),
  col_types = readr::cols(param1 = readr::col_double(),
                          param2 = readr::col_double()))
gi_dist <- distributional::dist_lognormal(
  mu = mean(gi_distribution_data$param1),
  sigma = mean(gi_distribution_data$param2))
gi <- as_discrete_pmf(gi_dist)
saveRDS(gi, 'tests/test_distributions/gi.rds')

## onset_to_notification
oldnotif_delay <- readRDS(paste0(not_synced_folder, '/ECDF_delay_constant_PCR.rds'))
min_delay <- 0
max_delay <- 41
onset_to_notification <- as_discrete_pmf(
  oldnotif_delay, min_delay, max_delay)
saveRDS(onset_to_notification, 'tests/test_distributions/onset_to_notification.rds')

