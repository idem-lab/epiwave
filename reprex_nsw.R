set.seed(1)

## real data
download.file("https://github.com/idem-lab/epiwave/raw/refs/heads/flat_prior/tests/test_distributions/cases_nsw_junoct.rds", case_tmp <- tempfile())
cases <- readRDS(case_tmp)
n_days <- length(cases)
days <- 1:n_days

# simulate case counts
car <- 0.75

download.file("https://github.com/idem-lab/epiwave/raw/refs/heads/main/tests/test_distributions/incubation.rds", inc_tmp <- tempfile())
incubation <- readRDS(inc_tmp)
download.file("https://github.com/idem-lab/epiwave/raw/refs/heads/main/tests/test_distributions/onset_to_notification.rds", ons_tmp <- tempfile())
onset_to_notification <- readRDS(ons_tmp)
inf_notif <- epiwave.params::add(incubation, onset_to_notification)

conv <- epiwave:::get_convolution_matrix(inf_notif, n = n_days)

# fit a deconvolution model with day of the week effect
library(greta)
incidence <- variable(0, dim = n_days)

expected_cases <- car * conv %*% incidence

# prior for dweek correction
dow_alpha <- greta::normal(1, 1,
                           truncation = c(0, Inf),
                           dim = c(1, 7))
dow_dist <- greta::dirichlet(dow_alpha,
                             n_realisations = 1)
dow_effect <- t(dow_dist * 7)

dow <- days %% 7 + 1
expected_cases_dow <- expected_cases * dow_effect[dow]
distribution(cases) <- poisson(expected_cases_dow)

m <- model(incidence)
draws <- mcmc(m)

# converges fine
rhats <- coda::gelman.diag(draws, autoburnin = FALSE, multivariate = FALSE)
max(rhats$psrf[, 1])

sry <- summary(draws)

n_sims <- 3
sims <- calculate(incidence, values = draws, nsim = n_sims)

ymax <- max(sry$quantiles[, "50%"])
plot(sry$statistics[, "Mean"] ~ days,
     ylab = "infection incidence",
     type = "n",
     ylim = c(0, ymax))
polygon(x = c(days, rev(days)),
        y = c(sry$quantiles[, "2.5%"],
              rev(sry$quantiles[, "97.5%"])),
        col = grey(0.9),
        lty = 0)
polygon(x = c(days, rev(days)),
        y = c(sry$quantiles[, "25%"],
              rev(sry$quantiles[, "75%"])),
        col = grey(0.7),
        lty = 0)
for (i in seq_len(n_sims)) {
  lines(sims$incidence[i, , 1],
        col = "red",
        lwd = 0.5)
}
lines(cases ~ days,
      lwd = 2)

dow_effect_est <- summary(calculate(dow_effect, values = draws))$statistics[, "Mean"]
# cbind(dow_effect_est, dow_effect_true)
