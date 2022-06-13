inv_pow_crt2 <- function(power, J, n, d_est = NULL, rho_est = NULL,
                         rsq2 = 0, K = 0, P = .5, alpha = .05,
                         test = "two.sided") {
  if (test == "two.sided") {
    if (is.null(J)) stop("Missing J")
    if (is.null(n)) stop("Missing n")
    df1 <- 1
    df2 <- J - K - 2
    cv <- stats::qf(1 - alpha, df1, df2)
    if (is.null(d_est)) {
      inv <- function(d_est) {
        ncp <- d_est^2 * (J * n * P * (1 - P) / (1 + (n * (1 - rsq2) - 1) * rho_est))
        stats::pf(cv, df1, df2, ncp, lower.tail = FALSE) - power
      }
      stats::uniroot(inv, c(0, 100))$root
    } else if (is.null(rho_est)) {
      inv <- function(rho_est) {
        ncp <- d_est^2 * (J * n * P * (1 - P) / (1 + (n * (1 - rsq2) - 1) * rho_est))
        stats::pf(cv, df1, df2, ncp, lower.tail = FALSE) - power
      }
      stats::uniroot(inv, c(0, 1))$root
    }
  } else if (test == "one.sided") {
    df <- J - K - 2
    cv <- stats::qt(1 - alpha, df)
    if (is.null(d_est)) {
      (cv - stats::qnorm(1 - power)) *
        sqrt((1 / (P * (1 - P)) * (1 + (n * (1 - rsq2) - 1) * rho_est)) / (J * n))
    } else if (is.null(rho_est)) {
      (J * n * P * (1 - P) * d_est^2 / (cv - stats::qnorm(1 - power))^2 - 1) /
        (n * (1 - rsq2) - 1)
    }
  }
}
