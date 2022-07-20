#' Assurance Level for Two-Level Multisite Randomized Trials
#'
#' \code{al_msrt2()} computes the assurance level over the specified uncertainty
#' about the parameters for a two-level MSRT design.
#'
#' @param d_est Effect size estimate, defined as
#'   \eqn{\delta = \frac{\gamma_{01}}{\tau^2 + \sigma^2}},
#'   where \eqn{\gamma_{01}} is the main effect of the treatment on the outcome,
#'   \eqn{\tau^2} is the variance of the cluster-specific random effect
#'   in the unconditional model (without covariates), and
#'   \eqn{\sigma^2} is the variance of the random error in the unconditional model.
#' @param d_sd Uncertainty level of the effect size estimate.
#' @param rho_est Intraclass correlation estimate, defined as
#'   \eqn{\rho = \frac{\tau^2}{\tau^2 + \sigma^2}}, where \eqn{\tau^2} and \eqn{\sigma^2}
#'   are the variance components in the unconditional model.
#' @param rho_sd Uncertainty level of the intraclass correlation estimate.
#' @param omega_est Estimate of the treatment effect hetereogeneity, defined as
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
#' @return The expected power given a two-level CRT design with J clusters each
#'   has n observations.
#' @export
#' @examples
#' al_crt2(J = 30, n = 100, d_est = .5, d_sd = .2, rho_est = .1, rho_sd = .05)

al_msrt2 <- function(J, n, d_est, d_sd, rho_est, rho_sd, omega_est, omega_sd,
                     rsq1 = 0, rsq2 = 0, K = 0, P = .5, power = .8, alpha = .05,
                     test = "two.sided") {
  d_est <- abs(d_est)

  # check if the multiplier of rho in the ncp equation is negative or positive
  # if negative, the higher the rho, the higher the ncp (lower.tail = FALSE)
  # if positive, the higher the rho, the lower the ncp (lower.tail = TRUE)
  if ((omega_est * (1 - rsq2) * P * (1 - P) * n + rsq1) < 1) {
    lower.tail <- FALSE
  } else {
    lower.tail <- TRUE
  }

  if (d_sd == 0) {
    if (rho_sd == 0) {
      if (omega_sd == 0) { # (1) d_sd = rho_sd = omega_sd = 0
        pow_msrt2(J = J, n = n, d_est = d_est, rho_est = rho_est,
                  omega_est = omega_est, rsq1 = rsq1, rsq2 = rsq2,
                  K = K, P = P, alpha = alpha, test = test)
      } else { # (2) d_sd = rho_sd = 0
        omega_ab <- gamma_ab(omega_est, omega_sd)
        stats::pgamma(
          inv_pow_msrt2(power = power, J = J, n = n, d_est = d_est,
                        rho_est = rho_est, rsq1 = rsq1, rsq2 = rsq2,
                        K = K, P = P, alpha = alpha, test = test),
          shape = omega_ab[1], rate = omega_ab[2]
        )
      }
    } else {
      if (omega_sd == 0) { # (3) d_sd = omega_sd = 0
        rho_ab <- get_ab(rho_est, rho_sd)
        stats::pbeta(
          inv_pow_msrt2(power = power, J = J, n = n, d_est = d_est,
                        omega_est = omega_est, rsq1 = rsq1, rsq2 = rsq2,
                        K = K, P = P, alpha = alpha, test = test),
          shape1 = rho_ab[1], shape2 = rho_ab[2], lower.tail = lower.tail
        )
      } else { # (4) d_sd = 0
        rho_ab <- get_ab(rho_est, rho_sd)
        omega_ab <- gamma_ab(omega_est, omega_sd)

        # boundary checking
        lo_bound <- inv_pow_msrt2(power = power, J = J, n = n, d_est = d_est,
                                  omega_est = 0, rsq1 = rsq1, rsq2 = rsq2,
                                  K = K, P = P, alpha = alpha, test = test)
        up_bound <- inv_pow_msrt2(power = power, J = J, n = n, d_est = d_est,
                                  omega_est = 1, rsq1 = rsq1, rsq2 = rsq2,
                                  K = K, P = P, alpha = alpha, test = test)
        if (lo_bound == 0 & up_bound == 0) lower.tail = FALSE

        cubature::cuhre(
          function(omega) {
            stats::pbeta(
              inv_pow_msrt2(power = power, J = J, n = n, d_est = d_est,
                            omega_est = omega, rsq1 = rsq1, rsq2 = rsq2,
                            K = K, P = P, alpha = alpha, test = test),
              shape1 = rho_ab[1], shape2 = rho_ab[2], lower.tail = lower.tail
            ) * stats::dgamma(omega, shape = omega_ab[1], rate = omega_ab[2])
          },
          lowerLimit = 0, upperLimit = 1
        )$integral
        # cubature::cuhre(
        #   function(rho) {
        #     stats::pgamma(
        #       inv_pow_msrt2(power = power, J = J, n = n, d_est = d_est,
        #                     rho_est = rho, rsq1 = rsq1, rsq2 = rsq2,
        #                     K = K, P = P, alpha = alpha, test = test),
        #       shape = omega_ab[1], rate = omega_ab[2], lower.tail = lower.tail
        #     ) * stats::dbeta(rho, shape1 = rho_ab[1], shape2 = rho_ab[2])
        #   },
        #   lowerLimit = 0, upperLimit = 1
        # )$integral
      }
    }
  } else {
    if (test == "two.sided") {
      if (rho_sd == 0) {
        if (omega_sd == 0) { # (5) rho_sd = omega_sd = 0
          d_star <- inv_pow_msrt2(power = power, J = J, n = n,
                                  rho_est = rho_est, omega_est = omega_est,
                                  rsq1 = rsq1, rsq2 = rsq2,
                                  K = K, P = P, alpha = alpha, test = test)
          stats::pnorm(d_star, mean = d_est, sd = d_sd, lower.tail = FALSE) +
            stats::pnorm(-d_star, mean = d_est, sd = d_sd, lower.tail = TRUE)
        } else { # (6) rho_sd = 0
          omega_ab <- gamma_ab(omega_est, omega_sd)
          cubature::cuhre(
            function(omega) {
              d_star <- inv_pow_msrt2(power = power, J = J, n = n,
                                      rho_est = rho_est, omega_est = omega,
                                      rsq1 = rsq1, rsq2 = rsq2,
                                      K = K, P = P, alpha = alpha, test = test)
              (stats::pnorm(d_star, mean = d_est, sd = d_sd, lower.tail = FALSE) +
                  stats::pnorm(-d_star, mean = d_est, sd = d_sd, lower.tail = TRUE)) *
                stats::dgamma(omega, shape = omega_ab[1], rate = omega_ab[2])
            },
            lowerLimit = 0, upperLimit = 1
          )$integral
        }
      } else {
        if (omega_sd == 0) { # (7) omega_sd = 0
          rho_ab <- get_ab(rho_est, rho_sd)
          cubature::hcubature(
            function(rho) {
              d_star <- inv_pow_msrt2(power = power, J = J, n = n,
                                      rho_est = rho, omega_est = omega_est,
                                      rsq1 = rsq1, rsq2 = rsq2,
                                      K = K, P = P, alpha = alpha, test = test)
              (stats::pnorm(d_star, mean = d_est, sd = d_sd, lower.tail = FALSE) +
                  stats::pnorm(-d_star, mean = d_est, sd = d_sd, lower.tail = TRUE)) *
                stats::dbeta(rho, shape1 = rho_ab[1], shape2 = rho_ab[2])
            },
            lowerLimit = 0, upperLimit = 1
          )$integral
        } else { # (8)
          rho_ab <- get_ab(rho_est, rho_sd)
          omega_ab <- gamma_ab(omega_est, omega_sd)
          cubature::cuhre(
            function(arg) {
              rho <- arg[1]
              omega <- arg[2]
              d_star <- inv_pow_msrt2(power = power, J = J, n = n,
                                      rho_est = rho, omega_est = omega,
                                      rsq1 = rsq1, rsq2 = rsq2,
                                      K = K, P = P, alpha = alpha, test = test)
              (stats::pnorm(d_star, mean = d_est, sd = d_sd, lower.tail = FALSE) +
                  stats::pnorm(-d_star, mean = d_est, sd = d_sd, lower.tail = TRUE)) *
                stats::dbeta(rho, shape1 = rho_ab[1], shape2 = rho_ab[2]) *
                stats::dgamma(omega, shape = omega_ab[1], rate = omega_ab[2])
            },
            lowerLimit = c(0, 0), upperLimit = c(1, 1)
          )$integral
        }
      }
    } else if (test == "one.sided") {
      if (rho_sd == 0) {
        if (omega_sd == 0) { # (5) rho_sd = omega_sd = 0
          d_star <- inv_pow_msrt2(power = power, J = J, n = n, d_est = d_est,
                                  rho_est = rho_est, omega_est = omega_est,
                                  rsq1 = rsq1, rsq2 = rsq2,
                                  K = K, P = P, alpha = alpha, test = test)
          stats::pnorm(d_star, mean = d_est, sd = d_sd, lower.tail = FALSE)
        } else { # (6) rho_sd = 0
          omega_ab <- gamma_ab(omega_est, omega_sd)
          cubature::cuhre(
            function(omega) {
              d_star <- inv_pow_msrt2(power = power, J = J, n = n, d_est = d_est,
                                      rho_est = rho_est, omega_est = omega,
                                      rsq1 = rsq1, rsq2 = rsq2,
                                      K = K, P = P, alpha = alpha, test = test)
              stats::pnorm(d_star, mean = d_est, sd = d_sd, lower.tail = FALSE) *
                stats::dgamma(omega, shape = omega_ab[1], rate = omega_ab[2])
            },
            lowerLimit = 0, upperLimit = 1
          )$integral
        }
      } else {
        if (omega_sd == 0) { # (7) omega_sd = 0
          rho_ab <- get_ab(rho_est, rho_sd)
          cubature::cuhre(
            function(rho) {
              d_star <- inv_pow_msrt2(power = power, J = J, n = n, d_est = d_est,
                                      rho_est = rho, omega_est = omega_est,
                                      rsq1 = rsq1, rsq2 = rsq2,
                                      K = K, P = P, alpha = alpha, test = test)
              stats::pnorm(d_star, mean = d_est, sd = d_sd, lower.tail = FALSE) *
                stats::dbeta(rho, shape1 = rho_ab[1], shape2 = rho_ab[2])
            },
            lowerLimit = 0, upperLimit = 1
          )$integral
        } else { # (8)
          rho_ab <- get_ab(rho_est, rho_sd)
          omega_ab <- gamma_ab(omega_est, omega_sd)
          cubature::cuhre(
            function(arg) {
              rho <- arg[1]
              omega <- arg[2]
              d_star <- inv_pow_msrt2(power = power, J = J, n = n, d_est = d_est,
                                      rho_est = rho, omega_est = omega,
                                      rsq1 = rsq1, rsq2 = rsq2,
                                      K = K, P = P, alpha = alpha, test = test)
              stats::pnorm(d_star, mean = d_est, sd = d_sd, lower.tail = FALSE) *
                stats::dbeta(rho, shape1 = rho_ab[1], shape2 = rho_ab[2]) *
                stats::dgamma(omega, shape = omega_ab[1], rate = omega_ab[2])
            },
            lowerLimit = 0, upperLimit = 1
          )$integral
        }
      }
    }
  }
}

