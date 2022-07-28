inv_se_msrt2 <- function(se, J, n, rho = NULL, omega = NULL,
                         rsq1 = 0, rsq2 = 0, K = 0, P = .5) {
  if (is.null(rho)) {
    inv <- function(rho_i) { # solve rho
      sqrt((rho_i * omega * (1 - rsq2) * P * (1 - P) * n +
              (1 - rho_i) * (1 - rsq1)) / (P * (1 - P) * J * n)) - se
    }
    # root finding & boundary checking
    inv_se_root(inv)
  } else if (is.null(omega)) {
    inv <- function(omega_i) {
      sqrt((rho * omega_i * (1 - rsq2) * P * (1 - P) * n +
              (1 - rho) * (1 - rsq1)) / (P * (1 - P) * J * n)) - se
    }
    # root finding & boundary checking
    inv_se_root(inv)
  }
}

inv_se_root <- function(inv, lb = 0, ub = 1) {
  root <- try(stats::uniroot(inv, c(lb, ub))$root,
              silent = TRUE)
  if (is(root, "try-error")) {
    if (inv(lb) > 0 & inv(ub) > 0) {
      # se > desired level for all icc/omega
      stop("Sample size is too small. ")
    } else if (inv(lb) < 0 & inv(ub) < 0) {
      1
    }
  } else {
    root
  }
}
