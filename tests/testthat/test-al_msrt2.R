test_that("Uncertainty in delta", {
  expect_equal(
    round(al_msrt2(J = 50, n = 30, d_est = .3, d_sd = .1,
                   rho_est = .2, rho_sd = 0,
                   omega_est = .3, omega_sd = 0), 4),
    .9115
  )
})

test_that("When there's uncertainty in rho", {
  expect_equal(
    round(al_msrt2(J = 20, n = 25, d_est = .3, d_sd = 0,
                   rho_est = .2, rho_sd = .1,
                   omega_est = .3, omega_sd = 0), 4),
    .8290
  )
  expect_equal(
    round(al_msrt2(J = 50, n = 7, d_est = .3, d_sd = 0,
                   rho_est = .2, rho_sd = .1,
                   omega_est = .3, omega_sd = 0), 4),
    .9716
  )
  expect_equal(
    al_msrt2(J = 50, n = 8, d_est = .3, d_sd = 0,
             rho_est = .2, rho_sd = .1,
             omega_est = .3, omega_sd = 0),
    1
  )
})

test_that("When there's uncertainty in omega", {
  expect_equal(
    round(al_msrt2(J = 20, n = 23, d_est = .3, d_sd = 0,
                   rho_est = .2, rho_sd = 0,
                   omega_est = .3, omega_sd = .1), 4),
    .5627
  )
})


test_that("When there's uncertainty in delta and rho", {
  expect_equal(
    round(al_msrt2(J = 50, n = 30, d_est = .3, d_sd = .1,
                   rho_est = .2, rho_sd = .1,
                   omega_est = .3, omega_sd = 0), 4),
    .9064
  )
})


test_that("When there's uncertainty in delta and omega", {
  expect_equal(
    round(al_msrt2(J = 50, n = 30, d_est = .3, d_sd = .1,
                   rho_est = .2, rho_sd = 0,
                   omega_est = .3, omega_sd = .1), 4),
    .9061
  )
})

test_that("When there's uncertainty in rho and omega", {
  expect_equal(
    round(al_msrt2(J = 20, n = 50, d_est = .3, d_sd = 0,
                   rho_est = .3, rho_sd = .1,
                   omega_est = .5, omega_sd = .1), 4),
    .449
  )
  expect_equal(
    round(al_msrt2(J = 100, n = 50, d_est = .3, d_sd = 0,
                   rho_est = .2, rho_sd = .2,
                   omega_est = .5, omega_sd = .1), 4),
    .9999
  )
})

test_that("When there's uncertainty in delta, rho, and omega", {
  expect_equal(
    round(al_msrt2(J = 20, n = 50, d_est = .3, d_sd = .1,
                   rho_est = .3, rho_sd = .1,
                   omega_est = .5, omega_sd = .1), 4),
    .4729
  )
})


