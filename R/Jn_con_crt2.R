Jn_con_crt2 <- function(d_est, rho_est, r2_est = 0, J = NULL, n = NULL,
                        K = 0, P = .5, power = .80, test = "two-tailed") {
  minJ <- K + 2 + 1
  lossJ <- function(J) {
    sum((pow_crt2(J = J, n = n, d_est = d_est, rho_est = rho_est,
                  r2_est = r2_est, test = test, P = P) - power)^2)
  }
  lossn <- function(n) {
    sum((pow_crt2(J = J, n = n, d_est = d_est, rho_est = rho_est,
                  r2_est = r2_est, test = test, P = P) - power)^2)
  }
  if (is.null(J)) {

    brent <- optim(K + 2 + 1, lossJ, lower = K + 2 + 1, upper = 1e6,
                   method = "Brent")
    # if L-BFGS-B does not converge, try using PORT routines
    if (brent$value > 1e-3) {
      port <- nlminb(minJ, lossJ, lower = minJ)
      # J <- optimize(lossJ, c(minJ, 1e6))$minimum
      if (port$objective > 1e-3) {
        lbfgsb <- optim(K + 2 + 1, lossJ, lower = K + 2 + 1, upper = Inf,
                        method = "L-BFGS-B")
        if (lbfgsb$value > 1e-3) {
          J <- lbfgsb$par
          warning(paste0("The algorithm fails to converge for the specified priors. ",
                         "There may not exist a solution for the desired expected ",
                         "power or assurance level. ",
                         "Please consider some lower power/assurance level. "))
        } else {
          J <- lbfgsb$par
        }
      } else {
        J <- port$par
      }
    } else {
      J <- brent$par
    }
    return(J)

  } else if (is.null(n)) {

    lbfgsb <- optim(1, lossn, lower = 1, upper = Inf, method = "L-BFGS-B")
    # if L-BFGS-B does not converge, try using PORT routines
    if (lbfgsb$value > 1e-3) {
      port <- nlminb(1, lossn, lower = 1)
      # if nlminb fails as well
      if (port$object > 1e-3) {
        stop(paste0("The algorithm fails to converge due to too few J ",
                    "for the specified priors. \n",
                    "Please consider raising J."))
      } else {
        n <- port$par
      }
    } else {
      n <- lbfgsb$par
    }
    return(n)
  }
}
