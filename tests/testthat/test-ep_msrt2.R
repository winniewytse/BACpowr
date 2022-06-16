test_that("Uncertainty in delta", {
  expect_equal(
    round(ep_msrt2(J = 50, n = 30, d_est = .3, d_sd = .1,
                   rho_est = .2, rho_sd = 0,
                   omega_est = .3, omega_sd = 0), 4),
    .9443
  )
})

test_that("When there's uncertainty in rho", {
  expect_equal(
    round(ep_msrt2(J = 20, n = 25, d_est = .3, d_sd = 0,
                   rho_est = .2, rho_sd = .1,
                   omega_est = .3, omega_sd = 0), 4),
    .8258
  )
  expect_equal(
    round(ep_msrt2(J = 50, n = 7, d_est = .3, d_sd = 0,
                   rho_est = .2, rho_sd = .1,
                   omega_est = .3, omega_sd = 0), 4),
    .8312
  )
  expect_equal(
    round(ep_msrt2(J = 50, n = 8, d_est = .3, d_sd = 0,
                   rho_est = .2, rho_sd = .1,
                   omega_est = .3, omega_sd = 0), 4),
    .8706
  )
})

test_that("When there's uncertainty in omega", {
  expect_equal(
    round(ep_msrt2(J = 20, n = 23, d_est = .3, d_sd = 0,
                   rho_est = .2, rho_sd = 0,
                   omega_est = .3, omega_sd = .1), 4),
    .8028
  )
})


test_that("When there's uncertainty in delta and rho", {
  expect_equal(
    round(ep_msrt2(J = 50, n = 30, d_est = .3, d_sd = .1,
                   rho_est = .2, rho_sd = .1,
                   omega_est = .3, omega_sd = 0), 4),
    .9413
  )
})


test_that("When there's uncertainty in delta and omega", {
  expect_equal(
    round(ep_msrt2(J = 50, n = 30, d_est = .3, d_sd = .1,
                   rho_est = .2, rho_sd = 0,
                   omega_est = .3, omega_sd = .1), 4),
    .9411
  )
})

test_that("When there's uncertainty in rho and omega", {
  expect_equal(
    round(ep_msrt2(J = 20, n = 50, d_est = .3, d_sd = 0,
                   rho_est = .3, rho_sd = .1,
                   omega_est = .5, omega_sd = .1), 4),
    .7793
  )
})

test_that("When there's uncertainty in delta, rho, and omega", {
  expect_equal(
    round(ep_msrt2(J = 20, n = 50, d_est = .3, d_sd = .1,
                   rho_est = .3, rho_sd = .1,
                   omega_est = .5, omega_sd = .1), 4),
    .7165
  )
})


