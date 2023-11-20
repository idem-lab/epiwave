estimate_reff <- function (infection_timeseries,
                           generation_interval_mass_fxns) {

    convolution_matrix <- get_convolution_matrix(
        generation_interval_mass_fxns, nrow(infection_timeseries))
    infectiousness <- convolution_matrix %*% infection_timeseries

    reff <- infection_timeseries / (infectiousness + .Machine$double.eps)

    greta_arrays <- module(
        infectiousness,
        reff)

    return(greta_arrays)
}
