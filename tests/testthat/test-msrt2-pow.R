#### Unit tests for the functions that help design two-level MSRTs ####

# pow_msrt2(), inv_pow_msrt2(), ep_msrt2(), al_msrt2(), Jn_msrt2()


# pow_msrt2() ------------------------------------------------------------------

test_that("Verying with PowerUpR", {
  expect_equal(
    round(pow_msrt2(J = 30, n = 10, delta = .3,
                    rho = .2, omega = 0), 4),
    .8016
    # PowerUpR::power.bira2(es = .3, rho2 = .2, omega2 = 0,
    #                       n = 10, J = 30)$power
  )
  expect_equal(
    round(pow_msrt2(J = 50, n = 20, delta = .2,
                    rho = .2, omega = .5), 4),
    .776
    # PowerUpR::power.bira2(es = .2, rho2 = .2, omega2 = .5,
    #                       n = 20, J = 50)$power
  )
})

# inv_pow_msrt2() --------------------------------------------------------------

test_that("Solve for delta (inverse power)", {
  expect_equal(round(
    inv_pow_msrt2(power = .8, J = 200, n = 50,
                  rho = .2, omega = .5, rsq2 = 0),
    4),
    0.0806)
})

test_that("Solve for rho (inverse power)", {
  expect_equal(
    round(inv_pow_msrt2(power = .8, J = 50, n = 20,
                        delta = .2, omega = .5, rsq2 = 0), 4),
    .1495
  )
  # checking
  # pow_msrt2(J = 50, n = 20, delta = .2, rho = .1495, omega = .5)
})

test_that("Solve for omega (inverse power)", {
  expect_equal(
    round(inv_pow_msrt2(power = .8, J = 30, n = 20,
                        delta = .3, rho = .2, rsq2 = 0), 4),
    .8066
  )
  # checking
  # pow_msrt2(J = 30, n = 20, delta = .3, rho = .2, omega = .8066)
})

test_that("Return ICC = 0 when power > desired level for all ICC (inverse power)", {
  expect_equal(inv_pow_msrt2(power = .8, J = 200, n = 50,
                             delta = .3, omega = .5, rsq2 = 0),
               0)
  # checking, power = 1
  # pow_msrt2(J = 200, n = 50, delta = .3, rho = 1, omega = .5)
})

test_that("Return omega = 0 when power > desired level for all omega (inverse power)", {
  expect_equal(inv_pow_msrt2(power = .8, J = 200, n = 50,
                             delta = .3, rho = .2, rsq2 = 0),
               0)
})

test_that("One-sided tests (inverse power)", {
  expect_equal(
    round(inv_pow_msrt2(power = .8, J = 200, n = 50, rho = .2,
                        omega = .5, rsq2 = 0, test = "one.sided"), 4),
    .0715
  )
  # pow_msrt2(J = 200, n = 50, delta = .0714, rho = .2,
  #           omega = .5, test = "one.sided")
  expect_equal(
    round(inv_pow_msrt2(power = .8, J = 50, n = 20, delta = .2,
                        omega = .5, rsq2 = 0, test = "one.sided"), 4),
    .3818
  )
  # checking
  # pow_msrt2(J = 50, n = 20, delta = .2, rho = .3819,
  #           omega = .5, test = "one.sided")
})


# ep_msrt2() -------------------------------------------------------------------

test_that("Uncertainty in delta (expected power)", {
  expect_equal(
    round(ep_msrt2(J = 50, n = 30, delta = .3, delta_sd = .1,
                   rho = .2, rho_sd = 0,
                   omega = .3, omega_sd = 0), 4),
    .9443
  )
})

test_that("Uncertainty in rho (expected power)", {
  expect_equal(
    round(ep_msrt2(J = 20, n = 25, delta = .3, delta_sd = 0,
                   rho = .2, rho_sd = .1,
                   omega = .3, omega_sd = 0), 4),
    .8258
  )
  expect_equal(
    round(ep_msrt2(J = 50, n = 7, delta = .3, delta_sd = 0,
                   rho = .2, rho_sd = .1,
                   omega = .3, omega_sd = 0), 4),
    .8312
  )
  expect_equal(
    round(ep_msrt2(J = 50, n = 8, delta = .3, delta_sd = 0,
                   rho = .2, rho_sd = .1,
                   omega = .3, omega_sd = 0), 4),
    .8706
  )
})

test_that("Uncertainty in omega (expected power)", {
  expect_equal(
    round(ep_msrt2(J = 20, n = 23, delta = .3, delta_sd = 0,
                   rho = .2, rho_sd = 0,
                   omega = .3, omega_sd = .1), 4),
    .8028
  )
})


test_that("Uncertainty in delta and rho (expected power)", {
  expect_equal(
    round(ep_msrt2(J = 50, n = 30, delta = .3, delta_sd = .1,
                   rho = .2, rho_sd = .1,
                   omega = .3, omega_sd = 0), 4),
    .9413
  )
})


test_that("Uncertainty in delta and omega (expected power)", {
  expect_equal(
    round(ep_msrt2(J = 50, n = 30, delta = .3, delta_sd = .1,
                   rho = .2, rho_sd = 0,
                   omega = .3, omega_sd = .1), 4),
    .9411
  )
})

test_that("Uncertainty in rho and omega (expected power)", {
  expect_equal(
    round(ep_msrt2(J = 20, n = 50, delta = .3, delta_sd = 0,
                   rho = .3, rho_sd = .1,
                   omega = .5, omega_sd = .1), 4),
    .7793
  )
})

test_that("Uncertainty in delta, rho, and omega (expected power)", {
  expect_equal(
    round(ep_msrt2(J = 20, n = 50, delta = .3, delta_sd = .1,
                   rho = .3, rho_sd = .1,
                   omega = .5, omega_sd = .1), 4),
    .7165
  )
})


# al_msrt2() -------------------------------------------------------------------

test_that("Uncertainty in delta (expected power)", {
  expect_equal(
    round(al_msrt2(J = 50, n = 30, delta = .3, delta_sd = .1,
                   rho = .2, rho_sd = 0,
                   omega = .3, omega_sd = 0), 4),
    .9115
  )
})

test_that("Uncertainty in rho (expected power)", {
  expect_equal(
    round(al_msrt2(J = 20, n = 25, delta = .3, delta_sd = 0,
                   rho = .2, rho_sd = .1,
                   omega = .3, omega_sd = 0), 4),
    .8290
  )
  expect_equal(
    round(al_msrt2(J = 50, n = 7, delta = .3, delta_sd = 0,
                   rho = .2, rho_sd = .1,
                   omega = .3, omega_sd = 0), 4),
    .9716
  )
  expect_equal(
    al_msrt2(J = 50, n = 8, delta = .3, delta_sd = 0,
             rho = .2, rho_sd = .1,
             omega = .3, omega_sd = 0),
    1
  )
})

test_that("Uncertainty in omega (expected power)", {
  expect_equal(
    round(al_msrt2(J = 20, n = 23, delta = .3, delta_sd = 0,
                   rho = .2, rho_sd = 0,
                   omega = .3, omega_sd = .1), 4),
    .5627
  )
})


test_that("Uncertainty in delta and rho (expected power)", {
  expect_equal(
    round(al_msrt2(J = 50, n = 30, delta = .3, delta_sd = .1,
                   rho = .2, rho_sd = .1,
                   omega = .3, omega_sd = 0), 4),
    .9064
  )
})


test_that("Uncertainty in delta and omega (expected power)", {
  expect_equal(
    round(al_msrt2(J = 50, n = 30, delta = .3, delta_sd = .1,
                   rho = .2, rho_sd = 0,
                   omega = .3, omega_sd = .1), 4),
    .9061
  )
})

test_that("Uncertainty in rho and omega (expected power)", {
  expect_equal(
    round(al_msrt2(J = 20, n = 50, delta = .3, delta_sd = 0,
                   rho = .3, rho_sd = .1,
                   omega = .5, omega_sd = .1), 4),
    .449
  )
  expect_equal(
    round(al_msrt2(J = 100, n = 50, delta = .3, delta_sd = 0,
                   rho = .2, rho_sd = .2,
                   omega = .5, omega_sd = .1), 4),
    .9999
  )
})

test_that("Uncertainty in delta, rho, and omega (expected power)", {
  expect_equal(
    round(al_msrt2(J = 20, n = 50, delta = .3, delta_sd = .1,
                   rho = .3, rho_sd = .1,
                   omega = .5, omega_sd = .1), 4),
    .4729
  )
})


# Jn_msrt2() -------------------------------------------------------------------

test_that("Without uncertainty", {
  expect_equal(
    Jn_msrt2(delta = .5, delta_sd = 0, rho = .1, rho_sd = 0,
             omega = .3, omega_sd = 0, n = 5)[1],
    26
    # PowerUpR::mrss.bira2(es = .5, rho2 = .1, omega2 = .3, n = 5)$J
  )
})

test_that("Uncertainty in omega", {
  expect_equal(
    Jn_msrt2(delta = .5, delta_sd = 0, rho = .1, rho_sd = 0,
             omega = .3, omega_sd = .1, n = 10)[1],
    15
  )
  # checking
  # ep_msrt2(J = 15, n = 10, delta = .5, delta_sd = 0, rho = .1, rho_sd = 0,
  #          omega = .3, omega_sd = .1)
})

test_that("Uncertainty in delta and rho", {
  expect_equal(
    Jn_msrt2(delta = .5, delta_sd = .1, rho = .1, rho_sd = .1,
             omega = .3, omega_sd = 0, n = 30)[1],
    8
  )
  # checking
  # ep_msrt2(J = 8, n = 30, delta = .5, delta_sd = .1, rho = .1, rho_sd = .1,
  #          omega = .3, omega_sd = 0)
  expect_equal(
    Jn_msrt2(delta = .5, delta_sd = .1, rho = .1, rho_sd = .1,
             omega = .3, omega_sd = 0, n = 30, al = .8)[1],
    10
  )
  # checking
  # al_msrt2(J = 10, n = 30, delta = .5, delta_sd = .1, rho = .1, rho_sd = .1,
  #          omega = .3, omega_sd = 0)
})

test_that("Uncertainty in rho and omega", {
  expect_equal(
    Jn_msrt2(delta = .5, delta_sd = 0, rho = .1, rho_sd = .1,
             omega = .3, omega_sd = .1, n = 5)[1],
    25
  )
  # checking
  # ep_msrt2(J = 25, n = 5, delta = .5, delta_sd = 0, rho = .1, rho_sd = .1,
  #          omega = .3, omega_sd = .1)
})

test_that("Uncertainty in delta and omega", {
  expect_equal(
    Jn_msrt2(delta = .5, delta_sd = .1, rho = .1, rho_sd = 0,
             omega = .3, omega_sd = .1, n = 5)[1],
    28
  )
  # checking
  # ep_msrt2(J = 28, n = 5, delta = .5, delta_sd = .1, rho = .1, rho_sd = 0,
  #          omega = .3, omega_sd = .1)
})

test_that("Uncertainty in delta, rho, and omega", {
  expect_equal(
    Jn_msrt2(delta = .5, delta_sd = .1, rho = .1, rho_sd = .1,
             omega = .3, omega_sd = .1, J = 30)[2],
    5
  )
  # checking
  # ep_msrt2(J = 30, n = 5, delta = .5, delta_sd = .1, rho = .1, rho_sd = .1,
  #          omega = .3, omega_sd = .1)
})




