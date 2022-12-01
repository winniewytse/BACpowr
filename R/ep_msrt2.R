#' Expected Power for Two-Level Multisite Randomized Trials
#'
#' \code{ep_msrt2()} computes the expected power over the specified uncertainty
#' about the parameters for a two-level MSRT design.
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
#' @param omega Estimate of the treatment effect hetereogeneity, defined as
#'   \eqn{\omega = \frac{\tau_1^2}{\tau_0^2}} where \eqn{\tau_0^2} is the variance of the
#'   intercept random component and \eqn{\tau_1^2} is the variance of the treatment
#'   random effect.
#' @param omega_sd Uncertainty level of the treatment effect hetereogeneity estimate.
#' @param rsq1 Estimate of variance explained by the level-1 (e.g., individual-level) covariates.
#' @param rsq2 Estimate of variance explained by the cluster-level covariates.
#' @param J Number of clusters. Determine \code{n} if \code{J} is specified.
#' @param n Cluster size. Determine \code{J} if \code{n} is specified.
#' @param K Number of cluster-level covariates.
#' @param P Proportion of the clusters that is treatment group.
#' @param power Desired statistical power to achieve. Default to be \code{.8}.
#' @param alpha Type I error rate. Default to be \code{.05}.
#' @param test One-sided or two-sided test. Options are either "one.sided" or "two.sided".
#' @param ... Additional arguments passed to \code{cubuture::hcubature()} to evaluate
#'   the integral.
#'
#' @return The expected power given a two-level CRT design with J clusters each
#'   with n observations.
#' @export
#' @examples
#' ep_msrt2(J = 30, n = 20, delta = .3, delta_sd = .1, rho = .2, rho_sd = .05,
#'          omega = .3, omega_sd = .05)

ep_msrt2 <- function(J, n, delta, delta_sd, rho, rho_sd, omega, omega_sd,
                     rsq1 = 0, rsq2 = 0, K = 0, P = .5, power = .8, alpha = .05,
                     test = "two.sided", ...) {

  # round extremely small delta_sd to 0 for computational stability
  if (delta_sd < .005) {delta_sd = 0} else {delta_sd = delta_sd}

  params <- list(J = J, n = n, rsq1 = rsq1, rsq2 = rsq2, K = K, P = P,
                 alpha = alpha, test = test)

  if (delta_sd == 0) {
    if (rho_sd == 0) {
      if (omega_sd == 0) { # (1) delta_sd = rho_sd = omega_sd = 0
        do.call(pow_msrt2,
                append(list(delta = delta, rho = rho, omega = omega), params))
      } else { # (2) delta_sd = rho_sd = 0
        omega_ab <- gamma_ab(omega, omega_sd)
        cubature::hcubature(
          function(x) {
            do.call(pow_msrt2,
                    append(list(delta = delta, rho = rho, omega = x), params)) *
              stats::dgamma(x, omega_ab[1], omega_ab[2])
          },
          lowerLimit = 0, upperLimit = 1,
          vectorInterface = TRUE, ...
        )$integral
      }
    } else {
      if (omega_sd == 0) { # (3) delta_sd = omega_sd = 0
        rho_ab <- beta_ab(rho, rho_sd)
        cubature::hcubature(
          function(x) {
            do.call(pow_msrt2,
                    append(list(delta = delta, rho = x, omega = omega), params)) *
              stats::dbeta(x, rho_ab[1], rho_ab[2])
          },
          lowerLimit = 0, upperLimit = 1,
          vectorInterface = TRUE, ...
        )$integral
      } else { # (4) delta_sd = 0
        rho_ab <- beta_ab(rho, rho_sd)
        omega_ab <- gamma_ab(omega, omega_sd)
        cubature::hcubature(
          function(matrix_arg) {
            matrix(apply(matrix_arg, 2, function(arg) {
              x <- arg[1]; y <- arg[2]
              do.call(pow_msrt2,
                      append(list(delta = delta, rho = x, omega = y), params)) *
                stats::dbeta(x, rho_ab[1], rho_ab[2]) *
                stats::dgamma(y, omega_ab[1], omega_ab[2])
            }), ncol = ncol(matrix_arg))
          },
          lowerLimit = c(0, 0), upperLimit = c(1, 1),
          vectorInterface = TRUE, ...
        )$integral
      }
    }
  } else {
    if (rho_sd == 0) {
      if (omega_sd == 0) { # (5) rho_sd = omega_sd = 0
        cubature::hcubature(
          function(x) {
            do.call(pow_msrt2,
                    append(list(delta = x, rho = rho, omega = omega), params)) *
              stats::dnorm(x, delta, delta_sd)
          },
          lowerLimit = -Inf, upperLimit = Inf,
          vectorInterface = TRUE, ...
        )$integral
      } else { # (6) rho_sd = 0
        omega_ab <- gamma_ab(omega, omega_sd)
        cubature::hcubature(
          function(matrix_arg) {
            matrix(apply(matrix_arg, 2, function(arg) {
              x <- arg[1]; y <- arg[2]
              do.call(pow_msrt2,
                      append(list(delta = x, rho = rho, omega = y), params)) *
                stats::dnorm(x, delta, delta_sd) *
                stats::dgamma(y, omega_ab[1], omega_ab[2])
            }), ncol = ncol(matrix_arg))
          },
          lowerLimit = c(-Inf, 0), upperLimit = c(Inf, 1),
          vectorInterface = TRUE, ...
        )$integral
      }
    } else {
      if (omega_sd == 0) { # (7) omega_sd = 0
        rho_ab <- beta_ab(rho, rho_sd)
        cubature::hcubature(
          function(matrix_arg) {
            matrix(apply(matrix_arg, 2, function(arg) {
              x <- arg[1]; y <- arg[2]
              do.call(pow_msrt2,
                      append(list(delta = x, rho = y, omega = omega), params)) *
                stats::dnorm(x, delta, delta_sd) *
                stats::dbeta(y, rho_ab[1], rho_ab[2])
            }), ncol = ncol(matrix_arg))
          },
          lowerLimit = c(-Inf, 0), upperLimit = c(Inf, 1),
          vectorInterface = TRUE, ...
        )$integral
      } else { # (8)
        rho_ab <- beta_ab(rho, rho_sd)
        omega_ab <- gamma_ab(omega, omega_sd)
        cubature::hcubature(
          function(matrix_arg) {
            matrix(apply(matrix_arg, 2, function(arg) {
              x <- arg[1]; y <- arg[2]; z <- arg[3]
              do.call(pow_msrt2,
                      append(list(delta = x, rho = y, omega = z), params)) *
                stats::dnorm(x, delta, delta_sd) *
                stats::dbeta(y, rho_ab[1], rho_ab[2]) *
                stats::dgamma(z, omega_ab[1], omega_ab[2])
            }), ncol = ncol(matrix_arg))
          },
          lowerLimit = c(-Inf, 0, 0), upperLimit = c(Inf, 1, 1),
          vectorInterface = TRUE, ...
        )$integral
      }
    }
  }
}
