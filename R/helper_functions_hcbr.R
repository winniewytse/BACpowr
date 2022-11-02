#' Define loss function for optimization
#' @param solve_for_J Indicate whether to define for the loss function for J.
#'   If FALSE, define it for n.
#' @param squared Indicate whether the loss function should be squared. If
#'   FALSE, set power_of to 1.
#' @export
define_loss <- function(solve_for_J = TRUE, squared = FALSE, ep, al,J, n, test,
                        reparameterize, list_params){

  power_of <- if (squared==TRUE) 2 else 1 #define what power to raise the loss

  library(zeallot) # where to put the library call? 'require zeallot'?

  c(d_est, d_sd, rho_est, rho_sd, rsq2, K, P,  power, alpha) %<-% list_params

  if(n == "") {n <- NULL}; if(J == "") {J <- NULL}

  # If AL is not specified, solve with EP (target) using ep_crt2().
  if (is.null(al) & !is.null(ep)) {
    criteria <- ep_crt2; target <- ep
  }
  # If EP is not specified, solve with AL (target) using al_crt2().
  if (is.null(ep) & !is.null(al)) {
    criteria <- al_crt2; target <- al
  }

  if(solve_for_J){
    loss <- function(J) {
      (criteria(J = J, n = n, d_est = d_est, d_sd = d_sd, rho_est = rho_est,
               rho_sd = rho_sd, rsq2 = rsq2, K = K, P = P, power = power,
               alpha = alpha, test = test,
               reparameterize = reparameterize) - target)^power_of
      }
    } else { # solve for n
      loss <- function(n) {
        (criteria(J = J, n = n, d_est = d_est, d_sd = d_sd, rho_est = rho_est,
                  rho_sd = rho_sd, rsq2 = rsq2, K = K, P = P, power = power,
                  alpha = alpha, test = test,
                  reparameterize = reparameterize) - target)^power_of
        }
    }
  return(loss)
  }


#' Compute intraclass correlation (ICC)
#'
#' \code{compute_icc()} computes ICC as
#'  \eqn{\rho = \frac{\tau^2}{\tau^2 + \sigma^2}}.
#'
#' @param sigma_sq The population within-group variance.
#' @param tau_sq The population between-group variance.
#' @return The ICC value, a number between 0 and 1.
#'
compute_icc <- function(r_sq, sigma_sq) {
  return(r_sq / (r_sq + sigma_sq))
}


# Solve Jn using the conventional approach
#' @export
Jn_crt2_c <- function(d_est, rho_est, rsq2 = 0,
                      J = NULL, n = NULL, K = 0, P = .5,
                      alpha = .05, power = .8, test = "two.sided",
                      reparameterize = FALSE) {

  if (is.null(J)) { # solve for J
    loss <- function(J) {
      pow_crt2(J = J, n = n, d_est = d_est, rho_est = rho_est,
               rsq2 = rsq2, test = test, P = P,
               reparameterize = reparameterize) - power
    }
    min <- K + 2 + 1
    J <- try(stats::uniroot(loss, c(min, 1e8))$root, silent = TRUE)
    if (class(J) == "try-error") {
      loss <- function(J) {
        (pow_crt2(J = J, n = n, d_est = d_est, rho_est = rho_est,
                  rsq2 = rsq2, test = test, P = P,
                  reparameterize = reparameterize) - power)^2
      }
      J <- Jn_optimize(start = min, loss = loss, lower = K + 3, upper = 1e6)
    }
  } else { # solve for n
    loss <- function(n) {
      pow_crt2(J = J, n = n, d_est = d_est, rho_est = rho_est,
               rsq2 = rsq2, test = test, P = P,
               reparameterize = reparameterize) - power
    }
    min <- 1
    n <- stats::uniroot(loss, c(min, 1e8))$root
  }
  return(cbind(J = J, n = n))
}


# Root finding & boundary checking
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

inv_pow_crt2 <- function(power, J, n, d_est = NULL, rho_est = NULL,
                         rsq2 = 0, K = 0, P = .5, alpha = .05,
                         test = "two.sided", reparameterize = FALSE) {
  df <- J - K - 2

  if (reparameterize) {
    # rho_est is defined as theta0 = tau^2 / sigma^2
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

inv_prec_msrt2 <- function(precision, J, n, rho = NULL, omega = NULL,
                           rsq1 = 0, rsq2 = 0, K = 0, P = .5, alpha = .05) {
  df <- J - K - 1
  t_crit <- qt(1 - alpha / 2, df = df)
  if (is.null(rho)) {
    inv <- function(rho_i) { # solve rho
      t_crit * sqrt((rho_i * omega * (1 - rsq2) * P * (1 - P) * n +
                       (1 - rho_i) * (1 - rsq1)) / (P * (1 - P) * J * n)) - precision
    }
    # root finding & boundary checking
    inv_prec_root(inv)
  } else if (is.null(omega)) {
    inv <- function(omega_i) {
      t_crit * sqrt((rho * omega_i * (1 - rsq2) * P * (1 - P) * n +
                       (1 - rho) * (1 - rsq1)) / (P * (1 - P) * J * n)) - precision
    }
    # root finding & boundary checking
    inv_prec_root(inv)
  }
}

inv_prec_root <- function(inv, lb = 0, ub = 1) {
  root <- try(stats::uniroot(inv, c(lb, ub))$root,
              silent = TRUE)
  if (is(root, "try-error")) {
    if (inv(lb) > 0 & inv(ub) > 0) {
      # precision > desired level for all icc/omega
      stop("Sample size is too small. ")
    } else if (inv(lb) < 0 & inv(ub) < 0) {
      1
    }
  } else {
    root
  }
}

