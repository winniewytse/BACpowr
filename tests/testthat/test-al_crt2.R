test_that("al_crt2() calculates the assurance level", {
  expect_equal(round(al_crt2(100, 3, .15, .39, .4, .14), 4),
               0.2284)
  expect_equal(round(al_crt2(19, 50, .8, .1, .2, .1), 4),
               0.7796)
  expect_error(Jn_crt2(.3, .1, .2, .1, J = 2))
  expect_equal(round(al_crt2(120, 20, .3, .1, .2, .2), 4),
               0.4985)
})
