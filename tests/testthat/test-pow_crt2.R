test_that("pow_crt2() calculate power of a two-level CRT", {
  expect_equal(round(pow_crt2(100, 50, .3, .2), 3),
               .892)
  expect_equal(round(pow_crt2(100, 50, .3, .2, P = .3), 3),
               .834)
  expect_lt(pow_crt2(100, 50, .3, .2, P = .8),
            pow_crt2(100, 50, .3, .2))
  expect_equal(pow_crt2(100, 50, .3, .2, P = .8),
               pow_crt2(100, 50, .3, .2, P = .2))
})

