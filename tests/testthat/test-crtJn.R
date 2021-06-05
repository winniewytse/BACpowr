pow1 <- crtJn(d_est = .5, d_sd = .2, rho_est = .1, rho_sd = .05, J = 30)
pow2 <- crtJn(d_est = .5, d_sd = .2, rho_est = .1, rho_sd = .05, n = 10)

test_that("crtJn() determines the minimum required J or n", {
  expect_equal(pow1,
               cbind(J = 30, n = 25))
})

test_that("crtJn() returns error with out of range input", {
  expect_error(crtJn(
    d_est = .5,
    d_sd = -.2,
    rho_est = .1,
    rho_sd = .05,
    J = 30
  ))
  expect_error(crtJn(
    d_est = .5,
    d_sd = .2,
    rho_est = -.1,
    rho_sd = .05,
    J = 30
  ))
  expect_error(crtJn(
    d_est = .5,
    d_sd = .2,
    rho_est = .1,
    rho_sd = -.05,
    J = 30
  ))
  expect_error(crtJn(
    d_est = .5,
    d_sd = .2,
    rho_est = .1,
    rho_sd = -.05,
    J = 0
  ))
})

# Based on the formula by McShane & Böckenholt
test_that("crtJn() gives results similar to single-level when ICC = 0", {
  expect_true(
    abs(crtJn(d_est = .5, d_sd = .2, rho_est = 0, rho_sd = 0, n = 1)[1] -
          182) <= 2
  )
})

test_that("crtJn() gives larger J/n with higher power", {
  expect_gt(
    crtJn(d_est = .5, d_sd = .2, rho_est = .1, rho_sd = .05, n = 10,
          power = .9)[1],
    pow2[1]
  )
  expect_gt(
    crtJn(d_est = .5, d_sd = .2, rho_est = .1, rho_sd = .05, J = 30,
          power = .81)[2],
    pow1[2]
  )
})

test_that("crtJn() gives larger J/n with two-tailed than one-tailed", {
  # To be added . . .
})
