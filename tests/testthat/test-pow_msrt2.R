

test_that("Verying with PowerUpR", {
  expect_equal(
    pow_msrt2(J = 30, n = 10, d_est = .3,
                    rho_est = 0.1966878, omega_est = 0),
    PowerUpR::power.bira2(es = .3, rho2 = 0.1966878, omega2 = 0,
                          n = 10, J = 30)$power
  )
  expect_equal(
    pow_msrt2(J = 50, n = 20, d_est = .2, rho_est = .1495, omega_est = .5),
    PowerUpR::power.bira2(es = .2, rho2 = .1495, omega2 = .5,
                          n = 20, J = 50)$power
  )
})


