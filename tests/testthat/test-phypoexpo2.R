# test if phypoexpo2 gives the same result as adding two exponentials
q <- rexp(1, 1)
rate_1 <- rexp(1, 1)
rate_2 <- rexp(1, 2)
nsim <- 1e7

test_that("hypoexponential distribution works", {
  expect_equal(mean(rexp(nsim, rate_1) + rexp(nsim, rate_2) <= q),
               phypoexpo2(q, rate_1, rate_2),
               tolerance = 1e-4)
})
