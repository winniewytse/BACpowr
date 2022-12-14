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


# Solve Jn using the conventional approach for a level-two CRT
#' @export
Jn_crt2_c <- function(delta, rho, rsq2 = 0, J = NULL, n = NULL, K = 0, P = .5,
                      alpha = .05, power = .8, test = "two.sided") {

  if (is.null(J)) { # solve for J
    loss <- function(J) {
      pow_crt2(J = J, n = n, delta = delta, rho = rho,
               rsq2 = rsq2, test = test, P = P) - power
    }
    min <- K + 2 + 1
    J <- try(stats::uniroot(loss, c(min, 1e8))$root, silent = TRUE)
    if (class(J) == "try-error") {
      loss <- function(J) {
        (pow_crt2(J = J, n = n, delta = delta, rho = rho,
                  rsq2 = rsq2, test = test, P = P) - power)^2
      }
      J <- optimize_Jn(start = min, loss = loss, lower = min, upper = 1e6,
                       solve = "J")
    }
  } else { # solve for n
    loss <- function(n) {
      pow_crt2(J = J, n = n, delta = delta, rho = rho, rsq2 = rsq2, test = test,
               P = P) - power
    }
    min <- 1
    n <- try(stats::uniroot(loss, c(min, 1e8))$root, silent = TRUE)
    if (class(n) == "try-error") {
      loss <- function(n) {
        (pow_crt2(J = J, n = n, delta = delta, rho = rho, rsq2 = rsq2,
                  test = test, P = P) - power)^2
      }
      n <- optimize_Jn(start = min, loss = loss, lower = min, upper = 1e6,
                       solve = "n")
    }
  }
  return(cbind(J = J, n = n))
}



# Solve Jn using the conventional approach for a level-two MSRT
#' @export
Jn_msrt2_c <- function(delta, rho, omega, rsq1 = 0, rsq2 = 0, J = NULL,
                       n = NULL, K = 0, P = .5, alpha = .05, power = .80,
                       test = "two.sided") {

  if (is.null(J)) { # solve J
    loss <- function(J) {
      pow_msrt2(J = J, n = n, delta = delta, rho = rho, omega = omega,
                rsq1 = rsq1, rsq2 = rsq2, test = test, K = K, P = P) - power
    }
    min <- K + 2 + 1
    J <- try(stats::uniroot(loss, c(min, 1e8))$root, silent = TRUE)
    if (class(J) == "try-error") {
      loss <- function(J) {
        (pow_msrt2(J = J, n = n, delta = delta, rho = rho, omega = omega,
                   rsq1 = rsq1, rsq2 = rsq2, test = test, K = K, P = P) - power)^2
      }
      J <- optimize_Jn(start = min, loss = loss, lower = K + 3, upper = 1e6)
    }
  } else { # solve n
    loss <- function(n) {
      pow_msrt2(J = J, n = n, delta = delta, rho = rho, omega = omega,
                rsq1 = rsq1, rsq2 = rsq2, test = test, K = K, P = P) - power
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
      1
    } else if(inv(lb) < 0 & inv(ub) < 0) {
      # power < desired level for all icc/omega
      stop("Effect size and/or sample size is too small. ")
    }
  } else {
    root
  }
}

# Inverse power function to solve delta or rho
inv_pow_crt2 <- function(power, J, n, delta = NULL, rho = NULL,
                         rsq2 = 0, K = 0, P = .5, alpha = .05,
                         test = "two.sided") {
  df <- J - K - 2

  # rho is defined as rho = tau^2 / (tau^2 + sigma^2)
  ncp <- function(delta, rho) {
    delta * sqrt(J * n * P * (1 - P) / (1 + (n * (1 - rsq2) - 1) * rho))
  }

  if (test == "two.sided") {
    cv <- stats::qt(1 - alpha / 2, df)
    if (is.null(delta)) {
      inv <- function(d) {
        stats::pt(cv, df, ncp = ncp(d, rho), lower.tail = FALSE) +
          stats::pt(-cv, df, ncp = ncp(d, rho), lower.tail = TRUE) - power
      }
      stats::uniroot(inv, c(0, 100))$root
    } else if (is.null(rho)) {
      inv <- function(rho) {
        stats::pt(cv, df, ncp = ncp(delta, rho), lower.tail = FALSE) +
          stats::pt(-cv, df, ncp = ncp(delta, rho), lower.tail = TRUE) - power
      }
      # root finding & boundary checking
      inv_pow_root(inv)
    }
  } else if (test == "one.sided") {
    cv <- stats::qt(1 - alpha, df)
    if (is.null(delta)) {
      inv <- function(d) {
        stats::pt(cv, df, ncp = ncp(d, rho), lower.tail = FALSE)  - power
      }
      stats::uniroot(inv, c(0, 100))$root
    } else if (is.null(rho)) {
      inv <- function(rho) {
        stats::pt(cv, df, ncp = ncp(delta, rho), lower.tail = FALSE) - power
      }
      # root finding & boundary checking
      inv_pow_root(inv)
    }
  }
}


inv_pow_msrt2 <- function(power, J, n, delta = NULL, rho = NULL, omega = NULL,
                          rsq1 = 0, rsq2 = 0, K = 0, P = .5, alpha = .05,
                          test = "two.sided") {

  # to avoid dividing 0
  # if (rho == 1 & omega == 0) rho <- .999

  df <- J - K - 1
  if (test == "two.sided") {
    cv <- stats::qt(1 - alpha / 2, df)
    if (is.null(delta)) {
      inv <- function(delta) { # solve d
        ncp <- delta * sqrt(P * (1 - P) * J * n /
                              (rho * omega * (1 - rsq2) * P * (1 - P) * n +
                                 (1 - rho) * (1 - rsq1)))
        stats::pt(cv, df, ncp, lower.tail = FALSE) +
          stats::pt(-cv, df, ncp, lower.tail = TRUE) - power
      }
      stats::uniroot(inv, c(0, 100))$root
    } else if (is.null(rho)) { # solve rho
      inv <- function(rho) {
        ncp <- delta * sqrt(P * (1 - P) * J * n /
                              (rho * omega * (1 - rsq2) * P * (1 - P) * n +
                                 (1 - rho) * (1 - rsq1)))
        stats::pt(cv, df, ncp, lower.tail = FALSE) +
          stats::pt(-cv, df, ncp, lower.tail = TRUE) - power
      }
      # root finding & boundary checking
      inv_pow_root(inv)
    } else if (is.null(omega)) { # solve omega
      inv <- function(omega) {
        ncp <- delta * sqrt(P * (1 - P) * J * n /
                              (rho * omega * (1 - rsq2) * P * (1 - P) * n +
                                 (1 - rho) * (1 - rsq1)))
        stats::pt(cv, df, ncp, lower.tail = FALSE) +
          stats::pt(-cv, df, ncp, lower.tail = TRUE) - power
      }
      # root finding & boundary checking
      inv_pow_root(inv)
    }
  } else if (test == "one.sided") {
    cv <- stats::qt(1 - alpha, df)
    if (is.null(delta)) {
      inv <- function(delta) { # solve d
        ncp <- delta * sqrt(P * (1 - P) * J * n /
                              (rho * omega * (1 - rsq2) * P * (1 - P) * n +
                                 (1 - rho) * (1 - rsq1)))
        stats::pt(cv, df, ncp, lower.tail = FALSE)  - power
      }
      stats::uniroot(inv, c(0, 100))$root
    } else if (is.null(rho)) { # solve rho
      inv <- function(rho) {
        ncp <- delta * sqrt(P * (1 - P) * J * n /
                              (rho * omega * (1 - rsq2) * P * (1 - P) * n +
                                 (1 - rho) * (1 - rsq1)))
        stats::pt(cv, df, ncp, lower.tail = FALSE) - power
      }
      # root finding & boundary checking
      inv_pow_root(inv)
    } else if (is.null(omega)) { # solve omega
      inv <- function(omega) {
        ncp <- delta * sqrt(P * (1 - P) * J * n /
                              (rho * omega * (1 - rsq2) * P * (1 - P) * n +
                                 (1 - rho) * (1 - rsq1)))
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

# Solve J/n by optimization method
optimize_Jn <- function(start, loss, lower, upper, message_par) {
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
        stop(do.call(err_message,
                     append(list(error = "optimize"), message_par)))
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

# check extreme cases
err_message <- function(given, goal, size, target, max_try = NULL, error) {
  if (given == "J") {
    m_given <- c("clusters", "participants", "number of clusters (J)")
  } else {
    m_given <- c("participants per cluster", "clusters",
                 "number of participants (n)")
  }
  if (goal == "ep") {
    m_goal <- c("expected power.", "expected power (`ep`)")
  } else {
    m_goal <- c("assurance level.", "assurace level (`al`)")
  }
  if (error == "max_try") {
    situation <- paste(
      "For a design with", size, m_given[1],
      "more than", max_try, m_given[2], "are needed to achieve",
      paste0(target * 100, "%"), m_goal[1]
    )
    additional <-
      "\nYou can also ask the program to try harder by increasing `max_try`."
  } else if (error == "optimize") {
    situation <- paste(
      "The algorithm fails to converge for a design with",
      size, paste0(m_given[1], ","),
      "meaning that there may not exist a solution for the number of",
      m_given[2], "that achieves", paste0(target * 100, "%"), m_goal[1]
    )
    additional <- NULL
  }
  suggestion <- paste(
    "Please consider increasing the", paste0(m_given[3], ","),
    "decreasing the desired", paste0(m_goal[2], ","),
    "or adjusting the level of uncertainty of each parameter."
  )
  paste0(situation, "\n", suggestion, additional)
}

Jn_try <- function(J, n, ep, al, params, max_try = 1e6, design) {
  if (design == "crt2") {
    ep_fun <- ep_crt2
    al_fun <- al_crt2
  } else if (design == "msrt2") {
    ep_fun <- ep_msrt2
    al_fun <- al_msrt2
  }
  if (!is.null(J)) { # solve n
    if (!is.null(ep)) {
      ep_try <- do.call(ep_fun, append(params, list(J = J, n = max_try)))
      if (ep_try < ep)
        stop(err_message("max_try", "J", "ep", J, ep, max_try))
    } else if (!is.null(al)) {
      al_try <- do.call(al_fun, append(params, list(J = J, n = max_try)))
      if (al_try < al)
        stop(err_message("max_try", "J", "al", J, al, max_try))
    }
  } else if (!is.null(n)) {
    if (!is.null(ep)) {
      ep_try <- do.call(ep_fun, append(params, list(J = max_try, n = n)))
      if (ep_try < ep)
        stop(err_message("max_try", "n", "ep", n, ep, max_try))
    } else if (!is.null(al)) {
      al_try <- do.call(al_fun, append(params, list(J = max_try, n = n)))
      if (al_try < al)
        stop(err_message("max_try", "n", "al", n, al, max_try))
    }
  }
}

