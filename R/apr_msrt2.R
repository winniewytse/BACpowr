#' Assurance Level of Precision for Two-Level Multisite Randomized Trials
#'
#' The probability of achieving the desired confidence interval half width or
#' narrower.
#'
#' @param rho Intraclass correlation value, defined as
#'   \eqn{\rho = \frac{\tau^2}{\tau^2 + \sigma^2}}, where \eqn{\tau^2} and \eqn{\sigma^2}
#'   are the variance components in the unconditional model.
#' @param rho_sd Uncertainty level of the intraclass correlation value.
#' @param omega Treatment effect hetereogeneity value, defined as
#'   \eqn{\omega = \frac{\tau_1^2}{\tau_0^2}} where \eqn{\tau_0^2} is the variance of the
#'   intercept random component and \eqn{\tau_1^2} is the variance of the treatment
#'   random effect.
#' @param omega_sd Uncertainty level of the treatment effect hetereogeneity value.
#' @param rsq1 Estimate of variance explained by the level-1 (e.g., individual-level) covariates.
#' @param rsq2 Estimate of variance explained by the cluster-level covariates.
#' @param J Number of clusters. Determine \code{n} if \code{J} is specified.
#' @param n Cluster size. Determine \code{J} if \code{n} is specified.
#' @param precision Desired confidence interval half width to achieve, defined as
#'        \eqn{t_{1 - \alpha / 2} \text{SE}(\hat \delta)}.
#' @param K Number of cluster-level covariates.
#' @param P Proportion of the clusters that is treatment group.
#' @return The expected confidence interval half width for two-level mulsite randomized trials.
#' @export
#'

apr_msrt2 <- function(rho, rho_sd, omega, omega_sd,
                      J, n, precision = 0.1, rsq1 = 0, rsq2 = 0,
                      K = 0, P = .5, alpha = .05, ...) {

  if (rho_sd == 0) {
    if (omega_sd == 0) { # (1) rho_sd = omega_sd = 0
      prec_msrt2(rho = rho, omega = omega, J = J, n = n,
                 rsq1 = rsq1, rsq2 = rsq2, K = K, P = P,
                 alpha = alpha)
    } else { # (2) rho_sd = 0
      omega_ab <- gamma_ab(omega, omega_sd)
      stats::pgamma(
        inv_prec_msrt2(precision = precision, J = J, n = n, rho = rho,
                       rsq1 = rsq1, rsq2 = rsq2, K = K, P = P,
                       alpha = alpha),
        shape = omega_ab[1], rate = omega_ab[2]
      )
    }
  } else {
    if (omega_sd == 0) { # (3) omega_sd = 0
      rho_ab <- get_ab(rho, rho_sd)
      stats::pbeta(
        inv_prec_msrt2(precision = precision, J = J, n = n, omega = omega,
                       rsq1 = rsq1, rsq2 = rsq2, K = K, P = P,
                       alpha = alpha),
        shape1 = rho_ab[1], shape2 = rho_ab[2]# , lower.tail = FALSE
      )
    } else { # (4)
      rho_ab <- get_ab(rho, rho_sd)
      omega_ab <- gamma_ab(omega, omega_sd)
      cubature::cuhre(
        function(x) {
          stats::pbeta(
            inv_prec_msrt2(precision = precision, J = J, n = n, omega = x,
                           rsq1 = rsq1, rsq2 = rsq2, K = K, P = P,
                           alpha = alpha),
            shape1 = rho_ab[1], shape2 = rho_ab[2]
          ) * stats::dgamma(x, shape = omega_ab[1], rate = omega_ab[2])
        },
        lowerLimit = 0, upperLimit = 1
      )$integral
    }
  }
}
