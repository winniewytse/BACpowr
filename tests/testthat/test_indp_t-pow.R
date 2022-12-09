
err_tolerance <- .0001

# pow_ind() --------------------------------------------------------------------

test_that("Verying with power.t.test()", {
  expect_lt(
    pow_indp_t(delta = .5, n1 = 50, n2 = 50) -
      power.t.test(50, delta = .5, type = "two.sample",
                   alternative = "two.sided")$pow,
    err_tolerance
  )
  expect_lt(
    pow_indp_t(delta = .2, n1 = 200, n2 = 200,
               test = "one.sided") -
      power.t.test(200, delta = .2, type = "two.sample",
                   alternative = "one.sided")$pow,
    err_tolerance
  )
})

# pow_inv2() -------------------------------------------------------------------

test_that("Inverse power function should return delta_L associated with L power", {
  expect_lt(
    pow_indp_t(
      delta = pow_inv2(power = .8, alpha = .05, n1 = 50, n2 = 50,
                       test = "two.sided"),
      n1 = 50, n2 = 50, test = "two.sided"
    ) - .8,
    err_tolerance
  )
  expect_lt(
    pow_indp_t(
      delta = pow_inv2(power = .6, alpha = .05, n1 = 50, n2 = 50,
                       test = "one.sided"),
      n1 = 50, n2 = 50, test = "one.sided"
    ) - .6,
    err_tolerance
  )
})

# ep_indp_t() ------------------------------------------------------------------

test_that("Expected power for independent sample t-test", {
  expect_lt(ep_indp_t(delta = .5, delta_sd = .2, n1 = 50, n2 = 50) - 0.6437375,
            err_tolerance)
  # sim_check_indp_t(n1 = 50, n2 = 50, delta = .5, delta_sd = .2)
})

# al_indp_t() ------------------------------------------------------------------

test_that("Assurance level for independent sample t-test", {
  expect_lt(al_indp_t(delta = .5, delta_sd = .2, n1 = 50, n2 = 50) - 0.370876,
            err_tolerance)
  # sim_check_indp_t(n1 = 50, n2 = 50, delta = .5, delta_sd = .2)
})

# n_indp_t() -------------------------------------------------------------------

test_that("Sample size for independent sample t-test", {
  expect_equal(n_indp_t(delta = .5, delta_sd = .2), 92)
  # ep_indp_t(delta = .5, delta_sd = .2, n1 = 92, n2 = 92)
  expect_equal(n_indp_t(delta = .5, delta_sd = .2, al = .8), 144)
  # ep_indp_t(delta = .5, delta_sd = .2, n1 = 92, n2 = 144)
})

