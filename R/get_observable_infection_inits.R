get_observable_infection_inits <- function (juris_no,
                                            obs_model_data) {

  idx_vals <- lapply(obs_model_data,
                     function(x) x$observable_idx_mat[, juris_no])
  inits_vals <- lapply(obs_model_data,
                       function(x) x$inits_values_mat[, juris_no])

  juris_specific_unique_idx <- sort(unique(
    as.vector(unlist(idx_vals))
  ))

  n_datasets <- length(idx_vals)

  juris_specific_inits <- vector(
    mode = 'numeric',
    length = length(juris_specific_unique_idx)
  )

  for (i in 1:length(juris_specific_inits)) {
    selector_idx <- juris_specific_unique_idx[i]

    inits_list <- list()

    for (dataset in 1:n_datasets) {

      dataset_inits <- get_inits(
        selector_idx,
        idx_vals[[dataset]],
        inits_vals[[dataset]])
      inits_list[dataset] <- dataset_inits

    }

    val <- mean(do.call(c, inits_list), na.rm = TRUE)

    juris_specific_inits[i] <- val

  }

  return(list(juris_specific_unique_idx = juris_specific_unique_idx,
              juris_specific_inits = juris_specific_inits))
}

get_inits <- function (selector_idx, observable_idx, inits_vals) {
  dataset_specific_selector_idx <- which(observable_idx == selector_idx)
  if (identical(dataset_specific_selector_idx, integer(0))) init <- NA
  else init <- inits_vals[dataset_specific_selector_idx]
  init
}
