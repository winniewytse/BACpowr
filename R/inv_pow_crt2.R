inv_pow_crt2 <- function(power, J, n, d_est = NULL, rho_est = NULL,
                         rsq2 = 0, K = 0, P = .5, alpha = .05,
                         test = "two.sided") {
  df <- J - K - 2
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
      # root finding & boundary checking
      inv_pow_root(inv)
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
      # root finding & boundary checking
      inv_pow_root(inv)
    }
  }
}

# root finding & boundary checking
inv_pow_root <- function(inv, lb = 0, ub = 1, tol = tol) {
  root <- try(stats::uniroot(inv, c(lb, ub), tol = tol)$root,
              silent = TRUE)
  if (is(root, "try-error")) {
    if (inv(lb) > 0 & inv(ub) > 0) {
      # power > desired level for all icc/omega
      0
    } else if(inv(lb) < 0 & inv(ub) < 0) {
      # power < desired level for all icc/omega
      stop("Effect size and/or sample size is too small. ")
    }
  } else {
    root
  }
}

