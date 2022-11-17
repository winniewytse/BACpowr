#' Expected Power for Two-Level CRTs
#'
#' \code{ep_crt2()} computes the expected power over the specified uncertainty
#' about the parameters for a two-level CRT design.
#'
#' @param delta Effect size estimate, defined as
#'   \eqn{\delta = \frac{\gamma_{01}}{\tau^2 + \sigma^2}},
#'   where \eqn{\gamma_{01}} is the main effect of the treatment on the outcome,
#'   \eqn{\tau^2} is the variance of the cluster-specific random effect
#'   in the unconditional model (without covariates), and
#'   \eqn{\sigma^2} is the variance of the random error in the unconditional model.
#' @param delta_sd Uncertainty level of the effect size estimate.
#' @param rho Intraclass correlation estimate, defined as
#'   \eqn{\rho = \frac{\tau^2}{\tau^2 + \sigma^2}}, where \eqn{\tau^2} and \eqn{\sigma^2}
#'   are the variance components in the unconditional model.
#' @param rho_sd Uncertainty level of the intraclass correlation estimate.
#' @param rsq2 Estimate of variance explained by the cluster-level covariates.
#' @param J Number of clusters. Determine \code{n} if \code{J} is specified.
#' @param n Cluster size. Determine \code{J} if \code{n} is specified.
#' @param K Number of cluster-level covariates.
#' @param P Proportion of the clusters that is treatment group.
#' @param power Desired statistical power to achieve. Default to be \code{.8}.
#' @param alpha Type I error rate. Default to be \code{.05}.
#' @param test One-sided or two-sided test. Options are either "one.sided" or "two.sided".
#' @param minEval Minimum number of evaluations of the integral.
#'   Defaults to `50` to stabilize results.
#' @param ... Additional arguments to be passed to \code{cubuture::cuhre()} to evaluate
#'   the integral.
#'
#' @return The expected power given a two-level CRT design with J clusters each
#'   with n observations.
#' @export
#' @examples
#' ep_crt2(J = 30, n = 100, delta = .5, delta_sd = .2, rho = .1, rho_sd = .05)
#' ep_crt2(J = 30, n = 100, delta = .5, delta_sd = .2, rho = .1, rho_sd = .05,
#'         rsq2 = .3)
ep_crt2 <- function(J, n, delta, delta_sd, rho, rho_sd,
                    rsq2 = 0, K = 0, P = .5, power = .8, alpha = .05,
                    test = "two.sided", minEval = 50, ...) {
  # round extremely small delta_sd to 0 for computational stability
  if (delta_sd < .005) {delta_sd = 0} else {delta_sd = delta_sd}

  if (delta_sd == 0) {
    if (rho_sd == 0) { # (1) delta_sd = rho_sd = 0
      pow_crt2(J = J, n = n, delta = delta, rho = rho,
               rsq2 = rsq2, K = K, P = P, alpha = alpha, test = test)
    } else {  # (2) delta_sd = 0
      rho_ab <- beta_ab(rho, rho_sd)
      cubature::cuhre(
        function(x) {
          pow_crt2(J = J, n = n, delta = delta, rho = x,
                   rsq2 = rsq2, K = K, P = P, alpha = alpha, test = test) *
            stats::dbeta(x, rho_ab[1], rho_ab[2])
        },
        lowerLimit = 0, upperLimit = 1,
        minEval = minEval, ...
      )$integral
    }
  } else {
    if (rho_sd == 0) { # (3) rho_sd = 0
      cubature::cuhre(
        function(x) {
          pow_crt2(J = J, n = n, delta = x, rho = rho,
                   rsq2 = rsq2, K = K, P = P, alpha = alpha, test = test) *
            stats::dnorm(delta, x, delta_sd)
        },
        lowerLimit = -Inf, upperLimit = Inf,
        minEval = minEval, ...
      )$integral
    } else {  # (4)
      rho_ab <- beta_ab(rho, rho_sd)
      cubature::cuhre(
        function(arg) {
          x <- arg[1]
          y <- arg[2]
          pow_crt2(J = J, n = n, delta = x, rho = y,
                   rsq2 = rsq2, K = K, P = P, alpha = alpha, test = test) *
            stats::dbeta(y, rho_ab[1], rho_ab[2]) *
            stats::dnorm(delta, x, delta_sd)
        },
        lowerLimit = c(-Inf, 0), upperLimit = c(Inf, 1),
        minEval = minEval, ...
      )$integral
    }
  }
}
