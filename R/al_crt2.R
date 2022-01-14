#' Assurance Level for Two-Level CRTs
#'
#' \code{ep_crt()} computes the assurance level of power given the estimates and
#' the uncertainty level of the parameter estimates for a two-level CRT.
#'
#' @param d_est Effect size estimate.
#' @param d_sd Uncertainty in the effect size estimate.
#' @param rho_est Intraclass correlation estimate.
#' @param rho_sd Uncertainty in the intraclass correlation estimate.
#' @param r2_est Estimate of variance explained by the cluster-level covariates.
#' @param r2_sd Uncertainty in the variance explained by the cluster-level covariates.
#' @param J Specified number of clusters.
#' @param n Specified cluster size.
#' @param K Number of cluster-level covariates.
#' @param power The desired level of power. Defaults to \code{.8}.
#' @param test One-tailed or two-tailed test.
#' @param abs.tol Absolute tolerance. Defaults to \code{1e-10}.
#' @param rel.tol Relative tolerance. Defaults to \code{1e-15}.
#' @return The expected power given certain J and n.
#' @export
#' @examples
#' al_crt2(J = 30, n = 100, d_est = .5, d_sd = .2, rho_est = .1, rho_sd = .05)
#' al_crt2(30, 100, .5, .2, .1, .05, .3, .1)

al_crt2 <- function(J, n, d_est, d_sd, rho_est, rho_sd,
                    r2_est = 0, r2_sd = 0, K = 0, power = .8,
                    test = "two-tailed",
                    abs.tol = 1e-50, rel.tol = 1e-3) {

  if (d_sd == 0) {
    if (rho_sd == 0) {
      if (r2_sd == 0) {             # (1) d_sd = rho_sd = r2_sd = 0
        pow_crt2(J, n, d_est, rho_est, r2_est, K, test)
      } else {                      # (2) d_sd = rho_sd = 0
        r2_ab <- get_ab(r2_est, r2_sd)
        1 - pbeta(inv_pow_crt2(power = power, J = J, n = n,
                               d_est = d_est, r2_est = r2_est),
                  shape1 = r2_ab[1], shape2 = r2_ab[2])
      }
    } else {
      if (r2_sd == 0) {             # (3) d_sd = r2_sd = 0
        rho_ab <- get_ab(rho_est, rho_sd)
        1 - pbeta(inv_pow_crt2(power = power, J = J, n = n,
                               d_est = d_est, r2_est = r2_est),
                  shape1 = rho_ab[1], shape2 = rho_ab[2])
      } else {                      # (4) d_sd = 0
        r2_ab <- get_ab(r2_est, r2_sd)
        rho_ab <- get_ab(rho_est, rho_sd)
        soluation <- cubature::cuhre(
          function(arg) {
            r2 <- arg[1]
            1 - pbeta(inv_pow_crt2(power = power, J = J, n = n,
                                   d_est = d_est, r2_est = r2_est),
                      shape1 = rho_ab[1], shape2 = rho_ab[2]) *
              stats::dbeta(r2, r2_ab[1], r2_ab[2])
          },
          lowerLimit = 0, upperLimit = 1,
          relTol = rel.tol, absTol = abs.tol
        )
      }
    }
  } else {
    if (rho_sd == 0) {
      if (r2_sd == 0) {             # (5) rho_sd = r2_sd = 0
        1 - pnorm(inv_pow_crt2(power = power, J = J, n = n,
                               rho_est = rho_est, r2_est = r2_est),
                  mean = d_est, sd = d_sd)
      } else {                      # (6) rho_sd = 0
        r2_ab <- get_ab(r2_est, r2_sd)
        solution <- cubature::cuhre(
          function(arg) {
            r2 <- arg[1]
            1 - pnorm(inv_pow_crt2(power = power, J = J, n = n,
                                   rho_est = rho_est, r2_est = r2_est),
                      mean = d_est, sd = d_sd) *
              stats::dbeta(r2, r2_ab[1], r2_ab[2])
          },
          lowerLimit = 0, upperLimit = 1,
          relTol = rel.tol, absTol = abs.tol
        )
      }
    } else {
      if (r2_sd == 0) {             # (7) r2_sd = 0
        rho_ab <- get_ab(rho_est, rho_sd)
        p <- 1 - pnorm(inv_pow_crt2(power = power, J = J, n = n,
                                    rho_est = rho_est, r2_est = r2_est),
                       mean = d_est, sd = d_sd)
        # when no power estimates are larger than the desired power level
        # assurance level = 0
        if (p == 0) return(0)
        solution <- cubature::cuhre(
          function(arg) {
            rho <- arg[1]
            p * stats::dbeta(rho, rho_ab[1], rho_ab[2])
          },
          lowerLimit = 0, upperLimit = 1,
          relTol = rel.tol, absTol = abs.tol
        )
      } else {                      # (8)
        r2_ab <- get_ab(r2_est, r2_sd)
        rho_ab <- get_ab(rho_est, rho_sd)
        solution <- cubature::cuhre(
          function(arg) {
            rho <- arg[1]
            r2 <- arg[2]
            1 - pnorm(inv_pow_crt2(power = power, J = J, n = n,
                                   rho_est = rho_est, r2_est = r2_est),
                      mean = d_est, sd = d_sd) *
              stats::dbeta(r2, r2_ab[1], r2_ab[2]) *
              stats::dbeta(rho, rho_ab[1], rho_ab[2])
          },
          lowerLimit = c(0, 0), upperLimit = c(1, 1),
          relTol = rel.tol, absTol = abs.tol
        )
      }
    }
  }
  # `prob` is the chisq probability that error is not a reliable estimate
  # as assurance level is too close to 0 given too small J or n
  if (solution$prob > 1.5e-2) return(0)
  solution$integral
}
