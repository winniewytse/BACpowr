test_that("Calculate the assurance level (two-sided)", {
  expect_equal(
    round(
      al_crt2(J = 100, n = 3, d_est = .15, d_sd = .39, rho_est = .4, rho_sd = .14),
      4),
    0.2939)
  expect_equal(
    round(
      al_crt2(J = 19, n = 50, d_est = .8, d_sd = .1, rho_est = .2, rho_sd = .1),
      4),
    0.7796)
  expect_equal(
    round(
      al_crt2(J = 120, n = 20, d_est = .3, d_sd = .1, rho_est = .2, rho_sd = .2),
      4),
    0.4985)
  expect_equal(
    round(
      al_crt2(J = 164, n = 20, d_est = .3, d_sd = 0, rho_est = .2, rho_sd = .18),
      4),
    0.7626)
})

test_that("Calculate the assurance level (one-sided)", {
  expect_equal(
    round(al_crt2(J = 100, n = 3, d_est = .15, d_sd = .39,
                  rho_est = .4, rho_sd = .14, test = "one.sided"), 4),
    .2693
  )
})
