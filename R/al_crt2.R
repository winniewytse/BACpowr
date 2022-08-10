#' Assurance Level for Two-Level CRTs
#'
#' \code{al_crt2()} computes the assurance level over the specified uncertainty
#' about the parameters for a two-level CRT design.
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
#' al_crt2(J = 30, n = 100, d_est = .5, d_sd = .2, rho_est = .1, rho_sd = .05,
#'         rsq2 = .3)

al_crt2 <- function(J, n, d_est, d_sd, rho_est, rho_sd,
                    rsq2 = 0, K = 0, P = .5, power = .8, alpha = .05,
                    test = "two.sided", reparameterize = FALSE) {
  d_est <- abs(d_est)

  if (rho_sd != 0) {
    if (reparameterize) {
      rho_p <- stats::pgamma
      rho_prior <- stats::dgamma
      rho_ab <- gamma_ab(rho_est, rho_sd)
      rho_up <- Inf
    } else {
      rho_p <- stats::pbeta
      rho_prior <- stats::dbeta
      rho_ab <- get_ab(rho_est, rho_sd)
      rho_up <- 1
    }
  }

  if (d_sd == 0) {
    if (rho_sd == 0) {              # (1) d_sd = rho_sd = 0
      pow_crt2(J = J, n = n, d_est = d_est, rho_est = rho_est,
               rsq2 = rsq2, K = K, P = P, alpha = alpha, test = test,
               reparameterize = reparameterize)
    } else {                        # (2) d_sd = 0
      # rho_ab <- get_ab(rho_est, rho_sd)
      rho_p(
        inv_pow_crt2(power = power, J = J, n = n,
                     d_est = d_est, rsq2 = rsq2, K = K, P = P,
                     alpha = alpha, test = test,
                     reparameterize = reparameterize),
        rho_ab[1], rho_ab[2]
      )
    }
  } else {
    if (test == "two.sided") {
      if (rho_sd == 0) {             # (3) rho_sd = 0
        d_star <- inv_pow_crt2(power = power, J = J, n = n,
                               rho_est = rho_est, rsq2 = rsq2, K = K, P = P,
                               alpha = alpha, test = test,
                               reparameterize = reparameterize)
        stats::pnorm(d_star, mean = d_est, sd = d_sd, lower.tail = FALSE) +
          stats::pnorm(-d_star, mean = d_est, sd = d_sd, lower.tail = TRUE)
      } else {                       # (4)
        # shapes <- get_ab(rho_est, rho_sd)
        cubature::cuhre(
          function(rho) {
            d_star <- inv_pow_crt2(power = power, J = J, n = n,
                                   rho_est = rho, rsq2 = rsq2, K = K, P = P,
                                   alpha = alpha, test = test,
                                   reparameterize = reparameterize)
            (stats::pnorm(d_star, mean = d_est, sd = d_sd, lower.tail = FALSE) +
                stats::pnorm(- d_star, mean = d_est, sd = d_sd, lower.tail = TRUE)) *
              rho_prior(rho, rho_ab[1], rho_ab[2])
          },
          lowerLimit = 0, upperLimit = rho_up
        )$integral
      }
    } else if (test == "one.sided") {
      if (rho_sd == 0) {             # (3) rho_sd = 0
        d_star <- inv_pow_crt2(power = power, J = J, n = n,
                               rho_est = rho_est, rsq2 = rsq2, K = K, P = P,
                               alpha = alpha, test = test,
                               reparameterize = reparameterize)
        stats::pnorm(d_star, mean = d_est, sd = d_sd, lower.tail = FALSE)
      } else {                       # (4)
        # shapes <- get_ab(rho_est, rho_sd)
        cubature::cuhre(
          function(rho) {
            d_star <- inv_pow_crt2(power = power, J = J, n = n,
                                   rho_est = rho, rsq2 = rsq2, K = K, P = P,
                                   alpha = alpha, test = test,
                                   reparameterize = reparameterize)
            stats::pnorm(d_star, mean = d_est, sd = d_sd, lower.tail = FALSE) *
              rho_prior(rho, rho_ab[1], rho_ab[2])
          },
          lowerLimit = 0, upperLimit = rho_up
        )$integral
      }
    }
  }
}
