# previous
n_dates <- 10
n_dates_observable <- 3
n_juris <- 2

infection_timeseries <- greta::zeros(n_dates, n_juris)
observable_idx <- matrix(c(0, 1, 1, 1, 0, 0), ncol = 2)
# observable_idx <- which(as.logical(observable_idx))

infection_observable <- greta::variable(lower = 0,
                                        dim = c(n_dates_observable,
                                                n_juris))

infection_timeseries[observable_idx,] <- infection_observable

# new
n_dates <- 10
n_juris <- 2

infection_timeseries <- greta::zeros(n_dates, n_juris)

observable_idx_mat <- matrix(c(0, 1, 1, 0, 0, 0,
                           0, 0, 0, 0, 1, 0,
                           0, 0, 0, 0, 0, 0,
                           0, 0), ncol = 2)
observable_idx <- apply(observable_idx_mat, 1:2, as.logical) #which(as.logical(observable_idx_mat))

no_true_observable_idx <- length(which(as.logical(observable_idx_mat)))

# flat_prior <- greta::variable(lower = 0, dim = c(n_dates, n_juris))
flat_prior <- greta::variable(lower = 0, dim = no_true_observable_idx)

# infection_timeseries[observable_idx] <- flat_prior[observable_idx]
infection_timeseries[observable_idx] <- flat_prior


