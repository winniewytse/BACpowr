
test_that("Solve omega that has SE = .05", {
  expect_equal(se_msrt2(rho = .3, J = 50, n = 30,
                        omega = inv_se_msrt2(se = .05, J = 50, n = 30,
                                             rho = .3)),
               .05)
})

test_that("Solve rho that has SE = .05", {
  expect_equal(se_msrt2(omega = .1, J = 50, n = 30,
                        rho = inv_se_msrt2(se = .05, J = 50, n = 30,
                                           omega = .1)),
               .05)
})
