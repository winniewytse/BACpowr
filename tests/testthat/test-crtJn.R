test_that("Jn_crt2() determines the minimum required J or n", {
  expect_equal(Jn_crt2(d_est = .5, d_sd = .2, rho_est = .1, rho_sd = .05, J = 30),
               cbind(J = 30, n = 25))
})
