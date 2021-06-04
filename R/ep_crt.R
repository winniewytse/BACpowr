#' Expected Power for Two-Level CRTs
#'
#' \code{ep_crt()} computes the expected power given the estimates of and the uncertainty
#' in the parameter estimates for a two-level CRT.
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
#' @param test One-tailed or two-tailed test.
#' @param abs.tol Absolute tolerance. Defaults to `1e-10`.
#' @param rel.tol Relative tolerance. Defaults to `1e-15`.
#' @return The expected power given certain J and n.
#' @export
#' @examples
#' ep_crt(J = 30, n = 100, d_est = .5, d_sd = .2, rho_est = .1, rho_sd = .05)
#' ep_crt(30, 100, .5, .2, .1, .05, .3, .1)
ep_crt <- function(J, n, d_est, d_sd, rho_est, rho_sd,
                   r2_est = 0, r2_sd = 0, K = 0,
                   test = "two-tailed",
                   abs.tol = 1e-50, rel.tol = 1e-3) {
  if (d_sd < .005) {d_sd = 0} else {d_sd = d_sd}
  if (d_sd == 0) {
    if (rho_sd == 0) {
      if (r2_sd == 0) {             # (1) d_sd = rho_sd = r2_sd = 0
        power_crt(J, n, d_est, rho_est, r2_est, K, test)
      } else {                      # (2) d_sd = rho_sd = 0
        r2_ab <- get_ab(r2_est, r2_sd)
        cubature::cuhre(
          function(arg) {
            r2 <- arg[1]
            power_crt(J, n, d_est, rho_est, r2, K, test) *
              stats::dbeta(r2, r2_ab[1], r2_ab[2])
          },
          lowerLimit = 0, upperLimit = 1,
          relTol = rel.tol, absTol= abs.tol
        )$integral
      }
    } else {
      if (r2_sd == 0) {             # (3) d_sd = r2_sd = 0
        rho_ab <- get_ab(rho_est, rho_sd)
        cubature::cuhre(
          function(arg) {
            rho <- arg[1]
            power_crt(J, n, d_est, rho, r2_est, K, test) *
              stats::dbeta(rho, rho_ab[1], rho_ab[2])
          },
          lowerLimit = 0, upperLimit = 1,
          relTol = rel.tol, absTol= abs.tol
        )$integral
      } else {                      # (4) d_sd = 0
        rho_ab <- get_ab(rho_est, rho_sd)
        r2_ab <- get_ab(r2_est, r2_sd)
        cubature::cuhre(
          function(arg) {
            rho <- arg[1]
            r2 <- arg[2]
            power_crt(J, n, d_est, rho, r2, test) *
              stats::dbeta(rho_ab[1], rho_ab[2]) *
              stats::dbeta(r2_ab[1], r2_ab[2])
          },
          lowerLimit = c(0, 0), upperLimit = c(1, 1),
          relTol = rel.tol, absTol= abs.tol
        )$integral
      }
    }
  } else {
    if (rho_sd == 0) {
      if (r2_sd == 0) {             # (5) rho_sd = r2_sd = 0
        cubature::cuhre(
          function(arg) {
            delta <- arg[1]
            power_crt(J, n, delta, rho_est, r2_est, K, test) *
              stats::dnorm(delta, d_est, d_sd)
          },
          lowerLimit = -Inf, upperLimit = Inf,
          relTol = rel.tol, absTol= abs.tol
        )$integral
      } else {                      # (6) rho_sd = 0
        r2_ab <- get_ab(r2_est, r2_sd)
        cubature::cuhre(
          function(arg) {
            delta <- arg[1]
            r2 <- arg[2]
            power_crt(J, n, delta, rho_est, r2, K, test) *
              stats::dbeta(r2, r2_ab[1], r2_ab[2]) *
              stats::dnorm(delta, d_est, d_sd)
          },
          lowerLimit = c(-Inf, 0), upperLimit = c(Inf, 1),
          relTol = rel.tol, absTol= abs.tol
        )$integral
      }
    } else {
      if (r2_sd == 0) {             # (7) r2_sd = 0
        rho_ab <- get_ab(rho_est, rho_sd)
        cubature::cuhre(
          function(arg) {
            delta <- arg[1]
            rho <- arg[2]
            power_crt(J, n, delta, rho, r2_est, K, test) *
              stats::dbeta(rho, rho_ab[1], rho_ab[2]) *
              stats::dnorm(delta, d_est, d_sd)
          },
          lowerLimit = c(-Inf, 0), upperLimit = c(Inf, 1),
          relTol = rel.tol, absTol= abs.tol
        )$integral
      } else {                      # (8)
        rho_ab <- get_ab(rho_est, rho_sd)
        r2_ab <- get_ab(r2_est, r2_sd)
        cubature::cuhre(
          function(arg) {
            delta <- arg[1]
            rho <- arg[2]
            r2 <- arg[3]
            power_crt(J, n, delta, rho, r2, K, test) *
              stats::dbeta(rho, rho_ab[1], rho_ab[2]) *
              stats::dbeta(r2, r2_ab[1], r2_ab[2]) *
              stats::dnorm(delta, d_est, d_sd)
          },
          lowerLimit = c(-Inf, 0, 0), upperLimit = c(Inf, 1, 1),
          relTol = rel.tol, absTol= abs.tol
        )$integral
      }
    }
  }
}
