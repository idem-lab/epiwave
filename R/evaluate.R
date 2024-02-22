# evaluate.timeseries_delay()
# evaluate.delay_probs()

evaluate.delay_probs <- function (day_diff, lookup, mass_val_col, delay_col) {

  lookup <- as.data.frame(lookup)
  day_diff[] <- lookup$massfun_combined[
    match(unlist(day_diff),
          lookup$total_delay)]
  day_diff[is.na(day_diff)] <- 0

  class(day_diff) <- 'numeric'
  day_diff

}

evaluate.timeseries_delay <- function () {

  con_list <- lapply(1:ncol(day_diff),
                            function (x) {
                              evaluate.delay_probs(
                                day_diff[,x], mass_functions[x,],
                                'massfun_combined', 'total_delay')
                            }
                     )
  con_mat <- do.call(cbind, con_list)
  con_mat

}


