test_that("Inverse power function", {
  expect_equal(round(inv_pow_crt2(.8, 200, 50, rho_est = .2, r2_est = 0), 4),
               0.1850)
  expect_equal(round(inv_pow_crt2(.8, 200, 50, d_est = .3, r2_est = 0), 4),
               0.5589)
})
