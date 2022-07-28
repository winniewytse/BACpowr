
#### se_msrt2 ####

test_that("Standard error for 2-level MSRT", {
  expect_lt(0.05773503 - se_msrt2(rho = .2, omega = .3, J = 50, n = 30),
            .0001)
})

#### inv_se_msrt2 ####

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

#### ese_msrt2 ####

test_that("Uncertainty in rho", {
  expect_lt(ese_msrt2(rho = .2, rho_sd = .1, omega = .3, omega_sd = 0, J = 50, n = 30) - 0.05867817,
            .0001)
})

test_that("Uncertainty in omega", {
  expect_lt(ese_msrt2(rho = .2, rho_sd = 0, omega = .3, omega_sd = .1, J = 50, n = 30) - 0.05867796,
            .0001)
})

test_that("Uncertainty in rho and omega", {
  expect_lt(ese_msrt2(rho = .2, rho_sd = .1, omega = .3, omega_sd = .1, J = 50, n = 30) - 0.059724,
            .0001)
})

#### ase_msrt2 ####

test_that("Uncertainty in rho", {
  expect_lt(ase_msrt2(rho = .2, rho_sd = .1, omega = .3, omega_sd = 0,
                      se = .06, J = 50, n = 30) - 0.6958234,
            .0001)
})

test_that("Uncertainty in omega", {
  expect_lt(ase_msrt2(rho = .2, rho_sd = 0, omega = .3, omega_sd = .1,
                      se = .06, J = 50, n = 30) - 0.6742127,
            .0001)
})

test_that("Uncertainty in rho and omega", {
  expect_lt(ase_msrt2(rho = .2, rho_sd = .1, omega = .3, omega_sd = .1,
                      se = .06, J = 50, n = 30) - 0.5960266,
            .0001)
})


#### Jn_msrt2 ####

test_that("Exepcted half width", {
  expect_equal(Jn_msrt2_se(rho = .2, rho_sd = .1, omega = .3, omega_sd = .1, J = 50, ese = .05)[2],
               62)
  # ese_msrt2(rho = .2, rho_sd = .1, omega = .3, omega_sd = .1, J = 50, n = 62)
  expect_equal(Jn_msrt2_se(rho = .2, rho_sd = .1, omega = .5, omega_sd = .3, n = 30, ese = .05)[1],
               72)
  # ese_msrt2(rho = .2, rho_sd = .1, omega = .5, omega_sd = .3, J = 72, n = 30)
})

test_that("Assurance level of half width", {
  expect_equal(Jn_msrt2_se(rho = .2, rho_sd = .1, omega = .5, omega_sd = .3, se = .06,
                           ase = .6, J = 50)[2],
               121)
  # ase_msrt2(rho = .2, rho_sd = .1, omega = .5, omega_sd = .3, se = .06, J = 50, n = 121)
  expect_equal(Jn_msrt2_se(rho = .2, rho_sd = .1, omega = .3, omega_sd = .1, n = 30,
                           se = .06, ase = .6)[1],
               51)
  # ase_msrt2(rho = .2, rho_sd = .1, omega = .3, omega_sd = .1, se = .06, J = 51, n = 30)
})
