test_that("Calculate the expected power", {
  # small delta_sd
  expect_equal(round(ep_crt2(J = 49, n = 20, delta = .4, delta_sd = .005,
                             rho = .2, rho_sd = 0, rsq2 = 0), 1),
               .8)
})

