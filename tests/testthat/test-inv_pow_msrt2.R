test_that("Solve for delta", {
  expect_equal(round(
    inv_pow_msrt2(power = .8, J = 200, n = 50,
                  rho_est = .2, omega_est = .5, rsq2 = 0),
    4),
    0.0806)
})

test_that("Solve for rho", {
  expect_equal(
    round(inv_pow_msrt2(power = .8, J = 50, n = 20,
                        d_est = .2, omega_est = .5, rsq2 = 0), 4),
    .1495
  )
  # checking
  # pow_msrt2(J = 50, n = 20, d_est = .2, rho_est = .1495, omega_est = .5)
})

test_that("Solve for omega", {
  expect_equal(
    round(inv_pow_msrt2(power = .8, J = 30, n = 20,
                        d_est = .3, rho_est = .2, rsq2 = 0), 4),
    .8066
  )
  # checking
  # pow_msrt2(J = 30, n = 20, d_est = .3, rho_est = .2, omega_est = .8066)
})

test_that("Return ICC = 1 when power > desired level for all ICC", {
  expect_equal(inv_pow_msrt2(power = .8, J = 200, n = 50,
                             d_est = .3, omega_est = .5, rsq2 = 0),
               1)
  # checking, power = 1
  # pow_msrt2(J = 200, n = 50, d_est = .3, rho_est = 1, omega_est = .5)
})

test_that("Return omega = 1 when power > desired level for all omega", {
  expect_equal(inv_pow_msrt2(power = .8, J = 200, n = 50,
                             d_est = .3, rho_est = .2, rsq2 = 0),
               1)
})

test_that("One-sided tests", {
  expect_equal(
    round(inv_pow_msrt2(power = .8, J = 200, n = 50, rho_est = .2,
                        omega_est = .5, rsq2 = 0, test = "one.sided"), 4),
    .0714
  )
  # pow_msrt2(J = 200, n = 50, d_est = .0714, rho_est = .2,
  #           omega_est = .5, test = "one.sided")
  expect_equal(
    round(inv_pow_msrt2(power = .8, J = 50, n = 20, d_est = .2,
                        omega_est = .5, rsq2 = 0, test = "one.sided"), 4),
    .3819
  )
  # checking
  # pow_msrt2(J = 50, n = 20, d_est = .2, rho_est = .3819,
  #           omega_est = .5, test = "one.sided")
})
