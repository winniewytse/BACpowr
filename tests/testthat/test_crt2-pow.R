#### Unit tests for the functions that help design two-level CRTs ####

# pow_crt2(), inv_pow_crt2(), ep_crt2(), al_crt2(), Jn_crt2()

# pow_crt2() -------------------------------------------------------------------

test_that("pow_crt2() calculate power of a two-level CRT", {
  expect_equal(round(pow_crt2(100, 50, .3, .2), 3),
               .892)
  expect_equal(round(pow_crt2(100, 50, .3, .2, P = .3), 3),
               .834)
  expect_lt(pow_crt2(100, 50, .3, .2, P = .8),
            pow_crt2(100, 50, .3, .2))
  expect_equal(pow_crt2(100, 50, .3, .2, P = .8),
               pow_crt2(100, 50, .3, .2, P = .2))
})


# inv_pow_crt2() ---------------------------------------------------------------

test_that("Calculate the effect size value when power is 80%", {
  expect_equal(round(
    inv_pow_crt2(power = .8, J = 200, n = 50,
                 rho = .2, rsq2 = 0),
    4),
    0.1850)
})

test_that("Calcualte the ICC value when power is 80%", {
  expect_equal(round(
    inv_pow_crt2(power = .8, J = 200, n = 50,
                 delta = .3, rsq2 = 0),
    4),
    0.5589)
})

test_that("Return ICC = 1 when power > desired level for all ICC", {
  expect_equal(inv_pow_crt2(power = .8, J = 200, n = 50,
                            delta = .5, rsq2 = 0),
               1)
  # checking, power = 1
  # pow_crt2(J = 200, n = 50, delta = .5, rho = 0, rsq2 = 0)
})

test_that("One-sided tests", {
  expect_equal(
    round(inv_pow_crt2(power = .8, J = 200, n = 50,
                       rho = .2, rsq2 = 0, test = "one.sided"), 4),
    .164
  )
  # checking
  # pow_crt2(J = 200, n = 50, delta = 0.164, rho = .2, test = "one.sided")
})

test_that("Special case", {
  expect_equal(
    inv_pow_crt2(power = .8, J = 1e6, n = 20, delta = .3, rho = NULL),
    1
  )
})


# ep_crt2() --------------------------------------------------------------------

test_that("Calculate the expected power", {
  # small delta_sd
  expect_equal(round(ep_crt2(J = 49, n = 20, delta = .4, delta_sd = .005,
                             rho = .2, rho_sd = 0, rsq2 = 0), 1),
               .8)
})

# al_crt2() --------------------------------------------------------------------

test_that("Calculate the assurance level (two-sided)", {
  expect_equal(
    round(
      al_crt2(J = 100, n = 3, delta = .15, delta_sd = .39, rho = .4, rho_sd = .14),
      4),
    0.2939)
  expect_equal(
    round(
      al_crt2(J = 19, n = 50, delta = .8, delta_sd = .1, rho = .2, rho_sd = .1),
      4),
    0.7796)
  expect_equal(
    round(
      al_crt2(J = 120, n = 20, delta = .3, delta_sd = .1, rho = .2, rho_sd = .2),
      4),
    0.4985)
  expect_equal(
    round(
      al_crt2(J = 164, n = 20, delta = .3, delta_sd = 0, rho = .2, rho_sd = .18),
      4),
    0.7626)
  expect_equal(al_crt2(1e6, 20, .3, 0, .2, .025), 1)
})

test_that("Calculate the assurance level (one-sided)", {
  expect_equal(
    round(al_crt2(J = 100, n = 3, delta = .15, delta_sd = .39,
                  rho = .4, rho_sd = .14, test = "one.sided"), 4),
    .2693
  )
})

# Jn_crt2() --------------------------------------------------------------------

test_that("Determine the required J for a desired expected power", {
  expect_equal(Jn_crt2(delta = .5, delta_sd = .2, rho = .1, rho_sd = .05,
                       J = 30, ep = .8),
               cbind(J = 30, n = 25))
})

test_that("Determine the reqruied J for a desired assurance level", {
  expect_lt(Jn_crt2(delta = .5, delta_sd = .1, rho = .2, rho_sd = .1,
                    rsq2 = 0, n = 50, al = .5)[1],
            Jn_crt2(delta = .25, delta_sd = .1, rho = .2, rho_sd = .1,
                    rsq2 = 0, n = 50, al = .5)[1])
  expect_equal(Jn_crt2(delta = .04436, delta_sd = .317,
                       rho = .01, rho_sd = .2294, n = 2, al = .8),
               cbind(J = 3132, n = 2))
})

test_that("Some special cases encountered before (assurance level)", {
  expect_equal(Jn_crt2(delta = .25, delta_sd = .1, rho = .2, rho_sd = .1,
                       rsq2 = 0, n = 50, al = .5),
               cbind(J = 123, n = 50))
  expect_equal(Jn_crt2(delta = .8, delta_sd = .1, rho = .2, rho_sd = .1,
                       rsq2 = 0, n = 50, al = .5)[1],
               15)
  expect_equal(Jn_crt2(delta = .3, delta_sd = .115, rho = .2, rho_sd = 0,
                       rsq2 = 0, n = 20, al = .8)[1],
               185)
  expect_equal(Jn_crt2(delta = .3, delta_sd = .23, rho = .2, rho_sd = 0,
                       n = 20, al = .8)[1],
               445)
})

test_that("Return an error when J isn't large enough", {
  expect_error(Jn_crt2(delta = .35, delta_sd = .1, rho = .15, rho_sd = .1,
                       J = 50))
})

test_that("Return an error when al doesn't achieve target level at max_try ", {
  expect_error(Jn_crt2(delta = .25, delta_sd = .1, rho = .2, rho_sd = .1,
                       rsq2 = 0, J = 50, al = .5))
})

test_that("Return an error if J smaller than the number of parameters", {
  expect_error(Jn_crt2(delta = .3, delta_sd = .1, rho = .2, rho_sd = .1,
                       J = 2, ep = .8))
})


test_that("When no uncertainty is specified", {
  expect_equal(Jn_crt2(delta = 0.958, delta_sd = 0, rho = 0.051, rho_sd = 0,
                       n = 3, ep = .8)[1],
               15)
  expect_equal(Jn_crt2(delta = 1.149, delta_sd = 0, rho = 0.258, rho_sd = 0,
                       n = 3, ep = .8)[1],
               15)
})

test_that("When only the uncertainty in rho is specified", {
  expect_equal(Jn_crt2(delta = .3, delta_sd = 0, rho = .2, rho_sd = .025,
                       n = 20, al = .8)[1],
               94)
})

test_that("With unequal treatment assignment", {
  expect_gt(Jn_crt2(delta = .3, delta_sd = .23, rho = .2, rho_sd = 0,
                    n = 20, P = .2, ep = .8)[1],
            Jn_crt2(delta = .3, delta_sd = .23, rho = .2, rho_sd = 0,
                    n = 20, ep = .8)[1])
  expect_gt(Jn_crt2(delta = .3, delta_sd = .23, rho = .2, rho_sd = 0,
                    n = 20, P = .2, al = .8)[1],
            Jn_crt2(delta = .3, delta_sd = .23, rho = .2, rho_sd = 0,
                    n = 20, al = .8)[1])
})


test_that("Solve Jn using the conventional approach", {
  expect_equal(round(Jn_crt2_c(delta = .5, rho = .2, n = 20)[1], 4),
               32.1782)
  expect_equal(Jn_crt2(delta = .5, delta_sd = 0, rho = .2, rho_sd = 0, n = 20)[1],
               33)
})


test_that("Change the desired statistical power", {
  expect_gt(Jn_crt2(delta = .5, delta_sd = .2, rho = .1, rho_sd = .05,
                    n = 30, power = .9, al = .6)[1],
            Jn_crt2(delta = .5, delta_sd = .2, rho = .1, rho_sd = .05,
                    n = 30, power = .8, al = .6)[1])
})

test_that("Print plots", {
  expect_is(Jn_crt2(delta = .5, delta_sd = .2, rho = .1, rho_sd = .05,
                    n = 30, power = .8, plot = TRUE)[[1]][[1]], "ggplot")
  expect_is(Jn_crt2(delta = .5, delta_sd = .2, rho = .1, rho_sd = .05,
                    n = 30, power = .8, al = .6, plot = TRUE)[[1]][[1]],
            "ggplot")
})

