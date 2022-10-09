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

# Solve J/n by optimization method
Jn_optimize <- function(start, loss, lower, upper, solve) {
  # try using PORT routines
  port <- stats::nlminb(start = start, objective = loss, lower = lower)
  # if PORT routines do not converge, try Brent
  if (port$objective > 1e-5) {
    brent <- try(stats::optim(par = start, fn = loss, lower = lower,
                              upper = upper, method = "Brent"),
                 silent = TRUE)
    # if Brent does not work, try LBFGSB
    if (class(brent) == "try-error") {
      condition <- "error"
    } else {
      if (brent$value > 1e-5) {
        condition <- "non-convergence"
      } else {
        condition <- "convergence"
      }
    }
    if (condition %in% c("error", "non-convergence")) {
      lbfgsb <- stats::optim(par = start, fn = loss, lower = lower,
                             upper = Inf, method = "L-BFGS-B")
      if (lbfgsb$value > 1e-5) {
        sol <- lbfgsb$par
        if (solve == "n") {
          warning("The algorithm fails to converge for the specified priors. ",
                  "Please consider increasing J or reducing the expected power
                  or", "assurance level. ")
        } else {
          warning(paste0("The algorithm fails to converge for the specified
                         priors. ", "There may not exist a solution for the
                         desired expected ", "power or assurance level. ",
                         "Please consider some lower power/assurance level."))
        }
      } else {
        sol <- lbfgsb$par
      }
    } else {
      sol <- brent$par
    }
  } else {
    sol <- port$par
  }
  return(sol)
}

# Solve Jn using the conventional approach
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
