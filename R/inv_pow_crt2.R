inv_pow_crt2 <- function(power, J, n, d_est = NULL, rho_est = NULL,
                         rsq2 = 0, K = 0, P = .5, alpha = .05,
                         test = "two.sided") {
  df <- J - K - 2
  ncp <- d_est * sqrt(J * n * P * (1 - P) / (1 + (n * (1 - rsq2) - 1) * rho_est))
  if (test == "two.sided") {
    cv <- stats::qt(1 - alpha / 2, df)
    if (is.null(d_est)) {
      inv <- function(d_est) {
        ncp <- d_est * sqrt(J * n * P * (1 - P) / (1 + (n * (1 - rsq2) - 1) * rho_est))
        stats::pt(cv, df, ncp, lower.tail = FALSE) +
          stats::pt(-cv, df, ncp, lower.tail = TRUE) - power
      }
      stats::uniroot(inv, c(0, 100))$root
    } else if (is.null(rho_est)) {
      inv <- function(rho_est) {
        ncp <- d_est * sqrt(J * n * P * (1 - P) / (1 + (n * (1 - rsq2) - 1) * rho_est))
        stats::pt(cv, df, ncp, lower.tail = FALSE) +
          stats::pt(-cv, df, ncp, lower.tail = TRUE) - power
      }
      stats::uniroot(inv, c(0, 1))$root
    }
  } else if (test == "one.sided") {
    cv <- stats::qt(1 - alpha, df)
    if (is.null(d_est)) {
      inv <- function(d_est) {
        ncp <- d_est * sqrt(J * n * P * (1 - P) / (1 + (n * (1 - rsq2) - 1) * rho_est))
        stats::pt(cv, df, ncp, lower.tail = FALSE)  - power
      }
      stats::uniroot(inv, c(0, 100))$root
    } else if (is.null(rho_est)) {
      inv <- function(rho_est) {
        ncp <- d_est * sqrt(J * n * P * (1 - P) / (1 + (n * (1 - rsq2) - 1) * rho_est))
        stats::pt(cv, df, ncp, lower.tail = FALSE) - power
      }
      stats::uniroot(inv, c(0, 1))$root
    }
  }
}
