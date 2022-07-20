#' Classical Statistical Power for Two-Level Multisite Randomized Trials
#'
#' \code{pow_msrt2()} computes the classical statistical power for a two-level MSRT.
#'
#' @param d_est Effect size estimate, defined as
#'   \eqn{\delta = \frac{\gamma_{01}}{\tau^2 + \sigma^2}},
#'   where \eqn{\gamma_{01}} is the main effect of the treatment on the outcome,
#'   \eqn{\tau^2} is the variance of the cluster-specific random effect
#'   in the unconditional model (without covariates), and
#'   \eqn{\sigma^2} is the variance of the random error in the unconditional model.
#' @param rho_est Intraclass correlation estimate, defined as
#'   \eqn{\rho = \frac{\tau^2}{\tau^2 + \sigma^2}}, where \eqn{\tau^2} and \eqn{\sigma^2}
#'   are the variance components in the unconditional model.
#' @param omega_est Estimate of the treatment effect hetereogeneity, defined as
#'   \eqn{\omega = \frac{\tau_1^2}{\tau_0^2}} where \eqn{\tau_0^2} is the variance of the
#'   intercept random component and \eqn{\tau_1^2} is the variance of the treatment
#'   random effect.
#' @param rsq1 Estimate of variance explained by the level-1 (e.g., individual-level) covariates.
#' @param rsq2 Estimate of variance explained by the level-2 (e.g., cluster-level) covariates.
#' @param J Number of clusters.
#' @param n Cluster size.
#' @param K Number of cluster-level covariates.
#' @param P Proportion of the sample that is treatment group.
#' @param alpha Type I error rate. Default to be \code{.05}.
#' @param test One-sided or two-sided test. Options are either "one.sided" or "two.sided".
#' @return Classical statistical power for a two-level CRT design with J clusters each has
#'   n observations.
#' @export
#' @examples
#' pow_msrt2(J = 30, n = 15, d_est = .3, rho_est = .2, omega_est = .3)
#' pow_msrt2(J = 10, n = 80, d_est = .3, rho_est = .2, omega_est = .2, P = .3)
pow_msrt2 <- function(J, n, d_est, rho_est, omega_est, rsq1 = 0, rsq2 = 0,
                      K = 0, P = .5, alpha = .05,
                      test = "two.sided") {
  if (J <= K + 1) stop(paste0("J needs to be larger than the number of parameters. ",
                              "Please increase J."))
  # if (rho_est == 1 & omega_est == 0)
  #   stop("The denominator of the ncp is 0 when rho = 1 and omega = 0.")

  df <- J - K - 1
  ncp <- d_est * sqrt(P * (1 - P) * J * n /
                        (rho_est * omega_est * (1 - rsq2) * P * (1 - P) * n +
                           (1 - rho_est) * (1 - rsq1)))
  if (test == "two.sided") {
    cv <- stats::qt(1 - alpha / 2, df)
    pow <- stats::pt(cv, df = df, ncp = ncp, lower.tail = FALSE) +
      stats::pt(-cv, df = df, ncp = ncp, lower.tail = TRUE)
  } else if (test == "one.sided") {
    cv <- stats::qt(1 - alpha, df)
    pow <- stats::pt(cv, df = df, ncp = ncp, lower.tail = FALSE)
  }
  return(pow)
}
