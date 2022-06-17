test_that("Without uncertainty", {
  expect_equal(
    Jn_msrt2(d_est = .5, d_sd = 0, rho_est = .1, rho_sd = 0,
             omega_est = .3, omega_sd = 0, n = 5)[1],
    26
    # PowerUpR::mrss.bira2(es = .5, rho2 = .1, omega2 = .3, n = 5)$J
  )
})

test_that("Uncertainty in omega", {
  expect_equal(
    Jn_msrt2(d_est = .5, d_sd = 0, rho_est = .1, rho_sd = 0,
             omega_est = .3, omega_sd = .1, n = 10)[1],
    15
  )
  # checking
  # ep_msrt2(J = 15, n = 10, d_est = .5, d_sd = 0, rho_est = .1, rho_sd = 0,
  #          omega_est = .3, omega_sd = .1)
})

test_that("Uncertainty in delta and rho", {
  expect_equal(
    Jn_msrt2(d_est = .5, d_sd = .1, rho_est = .1, rho_sd = .1,
             omega_est = .3, omega_sd = 0, n = 30)[1],
    8
  )
  # checking
  # ep_msrt2(J = 8, n = 30, d_est = .5, d_sd = .1, rho_est = .1, rho_sd = .1,
  #          omega_est = .3, omega_sd = 0)
  expect_equal(
    Jn_msrt2(d_est = .5, d_sd = .1, rho_est = .1, rho_sd = .1,
             omega_est = .3, omega_sd = 0, n = 30, al = .8)[1],
    10
  )
  # checking
  # al_msrt2(J = 10, n = 30, d_est = .5, d_sd = .1, rho_est = .1, rho_sd = .1,
  #          omega_est = .3, omega_sd = 0)
})

test_that("Uncertainty in rho and omega", {
  expect_equal(
    Jn_msrt2(d_est = .5, d_sd = 0, rho_est = .1, rho_sd = .1,
             omega_est = .3, omega_sd = .1, n = 5)[1],
    25
  )
  # checking
  # ep_msrt2(J = 25, n = 5, d_est = .5, d_sd = 0, rho_est = .1, rho_sd = .1,
  #          omega_est = .3, omega_sd = .1)
})

test_that("Uncertainty in delta and omega", {
  expect_equal(
    Jn_msrt2(d_est = .5, d_sd = .1, rho_est = .1, rho_sd = 0,
             omega_est = .3, omega_sd = .1, n = 5)[1],
    28
  )
  # checking
  # ep_msrt2(J = 28, n = 5, d_est = .5, d_sd = .1, rho_est = .1, rho_sd = 0,
  #          omega_est = .3, omega_sd = .1)
})

test_that("Uncertainty in delta, rho, and omega", {
  expect_equal(
    Jn_msrt2(d_est = .5, d_sd = .1, rho_est = .1, rho_sd = .1,
             omega_est = .3, omega_sd = .1, J = 30)[2],
    5
  )
  # checking
  # ep_msrt2(J = 30, n = 5, d_est = .5, d_sd = .1, rho_est = .1, rho_sd = .1,
  #          omega_est = .3, omega_sd = .1)
})

