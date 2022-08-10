test_that("Calculate the effect size value when power is 80%", {
  expect_equal(round(
    inv_pow_crt2(power = .8, J = 200, n = 50,
                 rho_est = .2, rsq2 = 0),
    4),
    0.1850)
})

test_that("Calcualte the ICC value when power is 80%", {
  expect_equal(round(
    inv_pow_crt2(power = .8, J = 200, n = 50,
                 d_est = .3, rsq2 = 0),
    4),
    0.5589)
})

test_that("Return ICC = 0 when power > desired level for all ICC", {
  expect_equal(inv_pow_crt2(power = .8, J = 200, n = 50,
                            d_est = .5, rsq2 = 0),
               0)
  # checking, power = 1
  # pow_crt2(J = 200, n = 50, d_est = .5, rho_est = 0, rsq2 = 0)
})

test_that("One-sided tests", {
  expect_equal(
    round(inv_pow_crt2(power = .8, J = 200, n = 50,
                       rho_est = .2, rsq2 = 0, test = "one.sided"), 4),
    .164
  )
  # checking
  # pow_crt2(J = 200, n = 50, d_est = 0.164, rho_est = .2, test = "one.sided")
})

