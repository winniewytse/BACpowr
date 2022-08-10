inv_pow_crt2 <- function(power, J, n, d_est = NULL, rho_est = NULL,
                         rsq2 = 0, K = 0, P = .5, alpha = .05,
                         test = "two.sided", reparameterize = FALSE) {
  df <- J - K - 2

  if (reparameterize) {
    # rho_est is defined as thata0 = tau^2 / sigma^2
    ncp <- function(d_est, rho_est) {
      d_est *
        sqrt(J * n * P * (1 - P) /
               (n * (1 - rsq2) * rho_est / (rho_est + 1) + (1 / (rho_est + 1))))
    }
  } else {
    # rho_est is defined as rho = tau^2 / (tau^2 + sigma^2)
    ncp <- function(d_est, rho_est) {
      d_est * sqrt(J * n * P * (1 - P) / (1 + (n * (1 - rsq2) - 1) * rho_est))
    }
  }

  if (test == "two.sided") {
    cv <- stats::qt(1 - alpha / 2, df)
    if (is.null(d_est)) {
      inv <- function(d) {
        # ncp <- d_est * sqrt(J * n * P * (1 - P) / (1 + (n * (1 - rsq2) - 1) * rho_est))
        stats::pt(cv, df, ncp = ncp(d, rho_est), lower.tail = FALSE) +
          stats::pt(-cv, df, ncp = ncp(d, rho_est), lower.tail = TRUE) - power
      }
      stats::uniroot(inv, c(0, 100))$root
    } else if (is.null(rho_est)) {
      inv <- function(rho) {
        # ncp <- d_est * sqrt(J * n * P * (1 - P) / (1 + (n * (1 - rsq2) - 1) * rho_est))
        stats::pt(cv, df, ncp = ncp(d_est, rho), lower.tail = FALSE) +
          stats::pt(-cv, df, ncp = ncp(d_est, rho), lower.tail = TRUE) - power
      }
      # root finding & boundary checking
      inv_pow_root(inv)
    }
  } else if (test == "one.sided") {
    cv <- stats::qt(1 - alpha, df)
    if (is.null(d_est)) {
      inv <- function(d) {
        # ncp <- d_est * sqrt(J * n * P * (1 - P) / (1 + (n * (1 - rsq2) - 1) * rho_est))
        stats::pt(cv, df, ncp = ncp(d, rho_est), lower.tail = FALSE)  - power
      }
      stats::uniroot(inv, c(0, 100))$root
    } else if (is.null(rho_est)) {
      inv <- function(rho) {
        # ncp <- d_est * sqrt(J * n * P * (1 - P) / (1 + (n * (1 - rsq2) - 1) * rho_est))
        stats::pt(cv, df, ncp = ncp(d_est, rho), lower.tail = FALSE) - power
      }
      # root finding & boundary checking
      inv_pow_root(inv)
    }
  }
}

# root finding & boundary checking
inv_pow_root <- function(inv, lb = 0, ub = 1) {
  root <- try(stats::uniroot(inv, c(lb, ub))$root,
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

