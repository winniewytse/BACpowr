test_that("Jn_crt2() determines the minimum required J or n", {
  expect_equal(Jn_crt2(d_est = .5, d_sd = .2, rho_est = .1, rho_sd = .05, J = 30),
               cbind(J = 30, n = 25))
  expect_equal(Jn_crt2(.25, .1, .2, .1, 0, 0, n = 50, al = .5),
               cbind(J = 111, n = 50))

  expect_equal(Jn_crt2(.8, .1, .2, .1, 0, 0, n = 50, al = .5)[1],
               236087)
  # very small d
  expect_warning(Jn_crt2(.05, .1, .2, .1, 0, 0, n = 3, al = .8))
})
