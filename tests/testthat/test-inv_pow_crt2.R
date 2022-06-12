test_that("Calculate the ICC value when power is 80%", {
  expect_equal(round(
    inv_pow_crt2(power = .8, J = 200, n = 50,
                 rho_est = .2, rsq2 = 0),
    4),
    0.1850)
})

test_that("Calcualte the effect size value when power is 80%", {
  expect_equal(round(
    inv_pow_crt2(power = .8, J = 200, n = 50,
                 d_est = .3, rsq2 = 0),
    4),
    0.5589)
})

