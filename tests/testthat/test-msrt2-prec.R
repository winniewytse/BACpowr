
#### prec_msrt2 ####

test_that("Standard error for 2-level MSRT", {
  expect_lt(0.1160229 - prec_msrt2(rho = .2, omega = .3, J = 50, n = 30),
            .0001)
})

#### inv_prec_msrt2 ####

test_that("Solve omega that has SE = .05", {
  expect_equal(prec_msrt2(rho = .3, J = 60, n = 25,
                          omega = inv_prec_msrt2(precision = .1,
                                                 J = 60, n = 25,
                                                 rho = .3)),
               .1)
})

test_that("Solve rho that has SE = .05", {
  expect_equal(prec_msrt2(omega = .1, J = 50, n = 30,
                          rho = inv_prec_msrt2(precision = .1,
                                               J = 50, n = 30,
                                               omega = .1)),
               .1)
})
#
# #### ese_msrt2 ####
#
# test_that("Uncertainty in rho", {
#   expect_lt(ese_msrt2(rho = .2, rho_sd = .1, omega = .3, omega_sd = 0, J = 50, n = 30) - 0.05867817,
#             .0001)
# })
#
# test_that("Uncertainty in omega", {
#   expect_lt(ese_msrt2(rho = .2, rho_sd = 0, omega = .3, omega_sd = .1, J = 50, n = 30) - 0.05867796,
#             .0001)
# })
#
# test_that("Uncertainty in rho and omega", {
#   expect_lt(ese_msrt2(rho = .2, rho_sd = .1, omega = .3, omega_sd = .1, J = 50, n = 30) - 0.059724,
#             .0001)
# })

#### apr_msrt2 ####

test_that("Uncertainty in rho", {
  expect_lt(apr_msrt2(rho = .2, rho_sd = .1, omega = .3, omega_sd = 0,
                      precision = .12, J = 50, n = 30) - 0.6630234,
            .0001)
})

test_that("Uncertainty in omega", {
  expect_lt(apr_msrt2(rho = .2, rho_sd = 0, omega = .3, omega_sd = .1,
                      precision = .12, J = 50, n = 30) - 0.6448033,
            .0001)
})

test_that("Uncertainty in rho and omega", {
  expect_lt(apr_msrt2(rho = .2, rho_sd = .1, omega = .3, omega_sd = .1,
                      precision = .12, J = 50, n = 30) - 0.5749551,
            .0001)
})


#### Jn_msrt2 ####

# test_that("Exepcted half width", {
#   expect_equal(Jn_msrt2_prec(rho = .2, rho_sd = .1, omega = .3,
#                              omega_sd = .1, J = 50, ese = .05)[2],
#                62)
#   # ese_msrt2(rho = .2, rho_sd = .1, omega = .3, omega_sd = .1, J = 50, n = 62)
#   expect_equal(Jn_msrt2_prec(rho = .2, rho_sd = .1, omega = .5,
#                              omega_sd = .3, n = 30, ese = .05)[1],
#                72)
#   # ese_msrt2(rho = .2, rho_sd = .1, omega = .5, omega_sd = .3, J = 72, n = 30)
# })

test_that("Assurance level of half width", {
  expect_equal(Jn_msrt2_prec(rho = .2, rho_sd = .1, omega = .5,
                             omega_sd = .3, precision = .12,
                             apr = .6, J = 50)[2],
               129)
  # apr_msrt2(rho = .2, rho_sd = .1, omega = .5, omega_sd = .3,
  #           precision = .12, J = 50, n = 129)
  expect_equal(Jn_msrt2_prec(rho = .2, rho_sd = .1, omega = .3,
                             omega_sd = .1, n = 30,
                             precision = .12, apr = .6)[1],
               51)
  # apr_msrt2(rho = .2, rho_sd = .1, omega = .3, omega_sd = .1,
  #           precision = .12, J = 51, n = 30)
})
