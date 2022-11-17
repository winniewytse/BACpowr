#' Assurance Level for Two-Level CRTs
#'
#' \code{al_crt2()} computes the assurance level over the specified uncertainty
#' about the parameters for a two-level CRT design.
#'
#' @param delta Effect size estimate, defined as
#'   \eqn{\delta = \frac{\gamma_{01}}{\tau^2 + \sigma^2}},
#'   where \eqn{\gamma_{01}} is the main effect of the treatment on the outcome,
#'   \eqn{\tau^2} is the variance of the cluster-specific random effect
#'   in the unconditional model (without covariates), and \eqn{\sigma^2} is the
#'   variance of the random error in the unconditional model.
#' @param delta_sd Uncertainty level of the effect size estimate.
#' @param rho Intraclass correlation estimate (ICC), defined as
#'   \eqn{\rho = \frac{\tau^2}{\tau^2 + \sigma^2}}, where \eqn{\tau^2} and
#'   \eqn{\sigma^2} are the variance components in the unconditional model.
#' @param rho_sd Uncertainty level of the intraclass correlation estimate.
#' @param rsq2 Estimate of variance explained by the cluster-level covariates.
#' @param J Number of clusters. Determine \code{n} if \code{J} is specified.
#' @param n Cluster size. Determine \code{J} if \code{n} is specified.
#' @param K Number of cluster-level covariates.
#' @param P Proportion of the clusters that belong to the treatment group.
#' @param power Desired statistical power. Defaults to \code{.8}.
#' @param alpha Type I error rate. Defaults to \code{.05}.
#' @param test One-sided or two-sided test. Defaults to "two.sided".
#' @return Assurance level given a two-level CRT design with J clusters with
#'    n observations.
#' @export
#' @examples
#' al_crt2(J = 30, n = 100, delta = .5, delta_sd = .2, rho = .1, rho_sd = .05)
#' al_crt2(J = 30, n = 100, delta = .5, delta_sd = .2, rho = .1, rho_sd = .05,
#'         rsq2 = .3)

al_crt2 <- function(J, n, delta, delta_sd, rho, rho_sd, rsq2 = 0, K = 0, P = .5,
                    power = .8, alpha = .05, test = "two.sided") {
  delta <- abs(delta)

  params <- list(J = J, n = n, rsq2 = rsq2, K = K, P = P, alpha = alpha,
                 test = test)

  if (delta_sd == 0) {
    if (rho_sd == 0) {              # (1) delta_sd = rho_sd = 0
      do.call(pow_crt2, append(list(rho = rho, delta = delta), params))
    } else {                        # (2) delta_sd = 0
      rho_ab <- beta_ab(rho, rho_sd)
      stats::pbeta(do.call(inv_pow_crt2,
                           append(list(power = power, delta = delta), params)),
                   rho_ab[1], rho_ab[2])
    }
  } else {
    if (test == "two.sided") {
      if (rho_sd == 0) {             # (3) rho_sd = 0
        d_L <- do.call(inv_pow_crt2, append(list(power = power, rho = rho), params))
        stats::pnorm(d_L, mean = delta, sd = delta_sd, lower.tail = FALSE) +
          stats::pnorm(-d_L, mean = delta, sd = delta_sd, lower.tail = TRUE)
      } else {                       # (4)
        rho_ab <- beta_ab(rho, rho_sd)
        cubature::cuhre(
          function(x) {
            d_L <- do.call(inv_pow_crt2, append(list(power = power, rho = x), params))
            (stats::pnorm(d_L, mean = delta, sd = delta_sd, lower.tail = FALSE) +
                stats::pnorm(- d_L, mean = delta, sd = delta_sd, lower.tail = TRUE)) *
              stats::dbeta(x, rho_ab[1], rho_ab[2])
          },
          lowerLimit = 0, upperLimit = 1
        )$integral
      }
    } else if (test == "one.sided") {
      if (rho_sd == 0) {             # (3) rho_sd = 0
        d_L <- do.call(inv_pow_crt2, append(list(power = power, rho = rho), params))
        stats::pnorm(d_L, mean = delta, sd = delta_sd, lower.tail = FALSE)
      } else {                       # (4)
        rho_ab <- beta_ab(rho, rho_sd)
        cubature::cuhre(
          function(x) {
            d_L <- do.call(inv_pow_crt2, append(list(power = power, rho = x), params))
            stats::pnorm(d_L, mean = delta, sd = delta_sd, lower.tail = FALSE) *
              stats::dbeta(x, rho_ab[1], rho_ab[2])
          },
          lowerLimit = 0, upperLimit = 1
        )$integral
      }
    }
  }
}
