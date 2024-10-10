## incubation period
incub_dist <- distributional::dist_weibull(
  shape = 1.83, scale = 4.93)
incubation <- parametric_dist_to_distribution(incub_dist)
saveRDS(incubation, 'tests/test_distributions/incubation.rds')

## generation interval distribution
not_synced_folder <- '../../BDSS_governance/BDSS/data'
gi_distribution_data <- readr::read_csv(
  file = paste0(not_synced_folder,'/nishiura_samples.csv'),
  col_types = readr::cols(param1 = readr::col_double(),
                          param2 = readr::col_double()))
gi_dist <- distributional::dist_lognormal(
  mu = mean(gi_distribution_data$param1),
  sigma = mean(gi_distribution_data$param2))
gi <- parametric_dist_to_distribution(gi_dist)
saveRDS(gi, 'tests/test_distributions/gi.rds')

## onset_to_notification
oldnotif_delay <- readRDS(paste0(not_synced_folder, '/ECDF_delay_constant_PCR.rds'))
min_delay <- 0
max_delay <- 41
onset_to_notification <- create_epiwave_massfun(
  min_delay, max_delay, oldnotif_delay, normalise = TRUE)
saveRDS(onset_to_notification, 'tests/test_distributions/onset_to_notification.rds')

