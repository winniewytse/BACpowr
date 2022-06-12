test_that("Calculate the expected power", {
  # small d_sd
  expect_equal(round(ep_crt2(J = 49, n = 20, d_est = .4, d_sd = .005,
                             rho_est = .2, rho_sd = 0, r2_est = 0,
                             r2_sd = 0), 1),
               .8)
})
