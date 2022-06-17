inv_pow_msrt2 <- function(power, J, n, d_est = NULL, rho_est = NULL, omega_est = NULL,
                          rsq1 = 0, rsq2 = 0, K = 0, P = .5, alpha = .05,
                          test = "two.sided") {

  # to avoid dividing 0
  # if (rho_est == 1 & omega_est == 0) rho_est <- .999

  df <- J - K - 1
  if (test == "two.sided") {
    cv <- stats::qt(1 - alpha / 2, df)
    if (is.null(d_est)) {
      inv <- function(d_est) { # solve d
        ncp <- d_est * sqrt(P * (1 - P) * J * n /
                              (rho_est * omega_est * (1 - rsq2) * P * (1 - P) * n +
                                 (1 - rho_est) * (1 - rsq1)))
        stats::pt(cv, df, ncp, lower.tail = FALSE) +
          stats::pt(-cv, df, ncp, lower.tail = TRUE) - power
      }
      stats::uniroot(inv, c(0, 100))$root
    } else if (is.null(rho_est)) { # solve rho
      inv <- function(rho_est) {
        ncp <- d_est * sqrt(P * (1 - P) * J * n /
                              (rho_est * omega_est * (1 - rsq2) * P * (1 - P) * n +
                                 (1 - rho_est) * (1 - rsq1)))
        stats::pt(cv, df, ncp, lower.tail = FALSE) +
          stats::pt(-cv, df, ncp, lower.tail = TRUE) - power
      }
      # root finding & boundary checking
      inv_pow_root(inv)
    } else if (is.null(omega_est)) { # solve omega
      inv <- function(omega_est) {
        ncp <- d_est * sqrt(P * (1 - P) * J * n /
                              (rho_est * omega_est * (1 - rsq2) * P * (1 - P) * n +
                                 (1 - rho_est) * (1 - rsq1)))
        stats::pt(cv, df, ncp, lower.tail = FALSE) +
          stats::pt(-cv, df, ncp, lower.tail = TRUE) - power
      }
      # root finding & boundary checking
      inv_pow_root(inv)
    }
  } else if (test == "one.sided") {
    cv <- stats::qt(1 - alpha, df)
    if (is.null(d_est)) {
      inv <- function(d_est) { # solve d
        ncp <- d_est * sqrt(P * (1 - P) * J * n /
                              (rho_est * omega_est * (1 - rsq2) * P * (1 - P) * n +
                                 (1 - rho_est) * (1 - rsq1)))
        stats::pt(cv, df, ncp, lower.tail = FALSE)  - power
      }
      stats::uniroot(inv, c(0, 100))$root
    } else if (is.null(rho_est)) { # solve rho
      inv <- function(rho_est) {
        ncp <- d_est * sqrt(P * (1 - P) * J * n /
                              (rho_est * omega_est * (1 - rsq2) * P * (1 - P) * n +
                                 (1 - rho_est) * (1 - rsq1)))
        stats::pt(cv, df, ncp, lower.tail = FALSE) - power
      }
      # root finding & boundary checking
      inv_pow_root(inv)
    } else if (is.null(omega_est)) { # solve omega
      inv <- function(omega_est) {
        ncp <- d_est * sqrt(P * (1 - P) * J * n /
                              (rho_est * omega_est * (1 - rsq2) * P * (1 - P) * n +
                                 (1 - rho_est) * (1 - rsq1)))
        stats::pt(cv, df, ncp, lower.tail = FALSE) - power
      }
      # root finding & boundary checking
      inv_pow_root(inv)
    }
  }
}

