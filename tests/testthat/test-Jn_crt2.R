test_that("Jn_crt2() determines the minimum required J or n", {
  expect_equal(Jn_crt2(d_est = .5, d_sd = .2, rho_est = .1, rho_sd = .05, J = 30),
               cbind(J = 30, n = 25))

  # comparison
  expect_lt(Jn_crt2(.5, .1, .2, .1, 0, 0, n = 50, al = .5)[1],
            Jn_crt2(.25, .1, .2, .1, 0, 0, n = 50, al = .5)[1])


  # special cases encountered before
  expect_equal(Jn_crt2(.25, .1, .2, .1, 0, 0, n = 50, al = .5),
               cbind(J = 123, n = 50))
  expect_equal(Jn_crt2(.8, .1, .2, .1, 0, 0, n = 50, al = .5)[1],
               15)
  expect_equal(Jn_crt2(.3, .115, .2, 0, n = 20, al = .8)[1],
               185)

  # very small d
  expect_warning(Jn_crt2(.05, .1, .2, .1, 0, 0, n = 3, al = .8))
  expect_equal(Jn_crt2(.05, .1, .2, .1, 0, 0, n = 3, al = .5)[1],
               6131)
  expect_error(Jn_crt2(.25, .1, .2, .1, 0, 0, J = 50, al = .5))

  # no uncertainty
  expect_equal(Jn_crt2(0.958, 0, 0.051, 0, n = 3)[1],
               15)
  expect_equal(Jn_crt2(1.149, 0, 0.258, 0, n = 3)[1],
               15)

  # no uncertainty in d
  expect_equal(Jn_crt2(.3, 0, .2, .025, n = 20, al = .8)[1],
               94)
})

