add_distributions <- function (delay1, delay2) {

  notif_delay_dist <- add(delay1, delay2)
  notif_full_delay_dist <- create_epiwave_massfun_timeseries(
    dates = infection_days,
    jurisdictions = jurisdictions,
    value = notif_delay_dist)
  notif_full_delay_dist

}
