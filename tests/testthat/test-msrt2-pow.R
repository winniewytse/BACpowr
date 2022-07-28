
#### pow_msrt2 ####

test_that("Verying with PowerUpR", {
  expect_equal(
    round(pow_msrt2(J = 30, n = 10, d_est = .3,
                    rho_est = .2, omega_est = 0), 4),
    .8016
    # PowerUpR::power.bira2(es = .3, rho2 = .2, omega2 = 0,
    #                       n = 10, J = 30)$power
  )
  expect_equal(
    round(pow_msrt2(J = 50, n = 20, d_est = .2,
                    rho_est = .2, omega_est = .5), 4),
    .776
    # PowerUpR::power.bira2(es = .2, rho2 = .2, omega2 = .5,
    #                       n = 20, J = 50)$power
  )
})

#### inv_pow_msrt2 ####

test_that("Solve for delta (inverse power)", {
  expect_equal(round(
    inv_pow_msrt2(power = .8, J = 200, n = 50,
                  rho_est = .2, omega_est = .5, rsq2 = 0),
    4),
    0.0806)
})

test_that("Solve for rho (inverse power)", {
  expect_equal(
    round(inv_pow_msrt2(power = .8, J = 50, n = 20,
                        d_est = .2, omega_est = .5, rsq2 = 0), 4),
    .1495
  )
  # checking
  # pow_msrt2(J = 50, n = 20, d_est = .2, rho_est = .1495, omega_est = .5)
})

test_that("Solve for omega (inverse power)", {
  expect_equal(
    round(inv_pow_msrt2(power = .8, J = 30, n = 20,
                        d_est = .3, rho_est = .2, rsq2 = 0), 4),
    .8066
  )
  # checking
  # pow_msrt2(J = 30, n = 20, d_est = .3, rho_est = .2, omega_est = .8066)
})

test_that("Return ICC = 0 when power > desired level for all ICC (inverse power)", {
  expect_equal(inv_pow_msrt2(power = .8, J = 200, n = 50,
                             d_est = .3, omega_est = .5, rsq2 = 0),
               0)
  # checking, power = 1
  # pow_msrt2(J = 200, n = 50, d_est = .3, rho_est = 1, omega_est = .5)
})

test_that("Return omega = 0 when power > desired level for all omega (inverse power)", {
  expect_equal(inv_pow_msrt2(power = .8, J = 200, n = 50,
                             d_est = .3, rho_est = .2, rsq2 = 0),
               0)
})

test_that("One-sided tests (inverse power)", {
  expect_equal(
    round(inv_pow_msrt2(power = .8, J = 200, n = 50, rho_est = .2,
                        omega_est = .5, rsq2 = 0, test = "one.sided"), 4),
    .0715
  )
  # pow_msrt2(J = 200, n = 50, d_est = .0714, rho_est = .2,
  #           omega_est = .5, test = "one.sided")
  expect_equal(
    round(inv_pow_msrt2(power = .8, J = 50, n = 20, d_est = .2,
                        omega_est = .5, rsq2 = 0, test = "one.sided"), 4),
    .3818
  )
  # checking
  # pow_msrt2(J = 50, n = 20, d_est = .2, rho_est = .3819,
  #           omega_est = .5, test = "one.sided")
})


#### ep_msrt2 ####

test_that("Uncertainty in delta (expected power)", {
  expect_equal(
    round(ep_msrt2(J = 50, n = 30, d_est = .3, d_sd = .1,
                   rho_est = .2, rho_sd = 0,
                   omega_est = .3, omega_sd = 0), 4),
    .9443
  )
})

test_that("Uncertainty in rho (expected power)", {
  expect_equal(
    round(ep_msrt2(J = 20, n = 25, d_est = .3, d_sd = 0,
                   rho_est = .2, rho_sd = .1,
                   omega_est = .3, omega_sd = 0), 4),
    .8258
  )
  expect_equal(
    round(ep_msrt2(J = 50, n = 7, d_est = .3, d_sd = 0,
                   rho_est = .2, rho_sd = .1,
                   omega_est = .3, omega_sd = 0), 4),
    .8312
  )
  expect_equal(
    round(ep_msrt2(J = 50, n = 8, d_est = .3, d_sd = 0,
                   rho_est = .2, rho_sd = .1,
                   omega_est = .3, omega_sd = 0), 4),
    .8706
  )
})

test_that("Uncertainty in omega (expected power)", {
  expect_equal(
    round(ep_msrt2(J = 20, n = 23, d_est = .3, d_sd = 0,
                   rho_est = .2, rho_sd = 0,
                   omega_est = .3, omega_sd = .1), 4),
    .8028
  )
})


test_that("Uncertainty in delta and rho (expected power)", {
  expect_equal(
    round(ep_msrt2(J = 50, n = 30, d_est = .3, d_sd = .1,
                   rho_est = .2, rho_sd = .1,
                   omega_est = .3, omega_sd = 0), 4),
    .9413
  )
})


test_that("Uncertainty in delta and omega (expected power)", {
  expect_equal(
    round(ep_msrt2(J = 50, n = 30, d_est = .3, d_sd = .1,
                   rho_est = .2, rho_sd = 0,
                   omega_est = .3, omega_sd = .1), 4),
    .9411
  )
})

test_that("Uncertainty in rho and omega (expected power)", {
  expect_equal(
    round(ep_msrt2(J = 20, n = 50, d_est = .3, d_sd = 0,
                   rho_est = .3, rho_sd = .1,
                   omega_est = .5, omega_sd = .1), 4),
    .7793
  )
})

test_that("Uncertainty in delta, rho, and omega (expected power)", {
  expect_equal(
    round(ep_msrt2(J = 20, n = 50, d_est = .3, d_sd = .1,
                   rho_est = .3, rho_sd = .1,
                   omega_est = .5, omega_sd = .1), 4),
    .7165
  )
})


#### al_msrt2 ####

test_that("Uncertainty in delta (expected power)", {
  expect_equal(
    round(al_msrt2(J = 50, n = 30, d_est = .3, d_sd = .1,
                   rho_est = .2, rho_sd = 0,
                   omega_est = .3, omega_sd = 0), 4),
    .9115
  )
})

test_that("Uncertainty in rho (expected power)", {
  expect_equal(
    round(al_msrt2(J = 20, n = 25, d_est = .3, d_sd = 0,
                   rho_est = .2, rho_sd = .1,
                   omega_est = .3, omega_sd = 0), 4),
    .8290
  )
  expect_equal(
    round(al_msrt2(J = 50, n = 7, d_est = .3, d_sd = 0,
                   rho_est = .2, rho_sd = .1,
                   omega_est = .3, omega_sd = 0), 4),
    .9716
  )
  expect_equal(
    al_msrt2(J = 50, n = 8, d_est = .3, d_sd = 0,
             rho_est = .2, rho_sd = .1,
             omega_est = .3, omega_sd = 0),
    1
  )
})

test_that("Uncertainty in omega (expected power)", {
  expect_equal(
    round(al_msrt2(J = 20, n = 23, d_est = .3, d_sd = 0,
                   rho_est = .2, rho_sd = 0,
                   omega_est = .3, omega_sd = .1), 4),
    .5627
  )
})


test_that("Uncertainty in delta and rho (expected power)", {
  expect_equal(
    round(al_msrt2(J = 50, n = 30, d_est = .3, d_sd = .1,
                   rho_est = .2, rho_sd = .1,
                   omega_est = .3, omega_sd = 0), 4),
    .9064
  )
})


test_that("Uncertainty in delta and omega (expected power)", {
  expect_equal(
    round(al_msrt2(J = 50, n = 30, d_est = .3, d_sd = .1,
                   rho_est = .2, rho_sd = 0,
                   omega_est = .3, omega_sd = .1), 4),
    .9061
  )
})

test_that("Uncertainty in rho and omega (expected power)", {
  expect_equal(
    round(al_msrt2(J = 20, n = 50, d_est = .3, d_sd = 0,
                   rho_est = .3, rho_sd = .1,
                   omega_est = .5, omega_sd = .1), 4),
    .449
  )
  expect_equal(
    round(al_msrt2(J = 100, n = 50, d_est = .3, d_sd = 0,
                   rho_est = .2, rho_sd = .2,
                   omega_est = .5, omega_sd = .1), 4),
    .9999
  )
})

test_that("Uncertainty in delta, rho, and omega (expected power)", {
  expect_equal(
    round(al_msrt2(J = 20, n = 50, d_est = .3, d_sd = .1,
                   rho_est = .3, rho_sd = .1,
                   omega_est = .5, omega_sd = .1), 4),
    .4729
  )
})





