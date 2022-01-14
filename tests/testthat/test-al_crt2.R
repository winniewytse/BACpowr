test_that("Jn_crt2() determines the minimum required J or n", {
  expect_lt(al_crt2(100, 3, .15, .39, .4, .14),
            0.2298606)
  expect_gt(al_crt2(19, 50, .8, .1, .2, .1),
            0.9515824)
})
