draws_msrt2 <- function(draws_d, draws_rho, draws_omega,
                        J, n, P = .5, K = 0, rsq1 = 0, rsq2 = 0,
                        test = "two.sided",
                        alpha = .05, power = .8, precision = .15,
                        goal = c("ep", "al", "apr")) {
  df <- J - 1
  cv <- qt(1 - alpha / 2, df)
  if (goal == "apr") {
    precision_draws <- cv * sqrt(
      (draws_rho * draws_omega * (1 - rsq2) *
         P * (1 - P) * n + (1 - draws_rho) * (1 - rsq1)) /
        (P * (1 - P) * J * n)
    )
    mean(precision_draws < precision)
  } else {
    draws_ncp <- draws_d * sqrt(
      P * (1 - P) * J * n /
        (draws_rho * draws_omega * (1 - rsq2) * P * (1 - P) * n +
           (1 - draws_rho) * (1 - rsq1))
    )
    draws_pow <- pt(cv, df = df, ncp = draws_ncp, lower.tail = FALSE) +
      pt(-cv, df = df, ncp = draws_ncp, lower.tail = TRUE)
    if (goal == "ep") {
      mean(draws_pow)
    } else if (goal == "al") {
      mean(draws_pow > power)
    }
  }
}
