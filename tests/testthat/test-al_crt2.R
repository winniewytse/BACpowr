test_that("al_crt2() calculates the assurance level", {
  expect_lt(al_crt2(100, 3, .15, .39, .4, .14),
            0.2298606)
  expect_gt(al_crt2(19, 50, .8, .1, .2, .1),
            0.9515824)
  expect_error(Jn_crt2(.3, .1, .2, .1, J = 2))
})
