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

test_that("Return warning when J isn't large enough", {
  expect_warning(Jn_crt2(delta = .35, delta_sd = .1, rho = .15, rho_sd = .1,
                         J = 50))
})

test_that("When the effect size is very small (assurance level)", {
  expect_equal(Jn_crt2(delta = .05, delta_sd = .1, rho = .2, rho_sd = .1,
                       rsq2 = 0, n = 3, al = .5)[1],
               2634)
  expect_error(Jn_crt2(delta = .25, delta_sd = .1, rho = .2, rho_sd = .1,
                       rsq2 = 0, J = 50, al = .5))
})

test_that("Return error if J is not large enough", {
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
