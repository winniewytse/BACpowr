#' Precision for Two-Level Multisite Randomized Trials
#'
#' Precision is defined as the confidence interval half width.
#'
#' @param rho Intraclass correlation value, defined as
#'   \eqn{\rho = \frac{\tau^2}{\tau^2 + \sigma^2}}, where \eqn{\tau^2} and \eqn{\sigma^2}
#'   are the variance components in the unconditional model.
#' @param omega Treatment effect hetereogeneity value, defined as
#'   \eqn{\omega = \frac{\tau_1^2}{\tau_0^2}} where \eqn{\tau_0^2} is the variance of the
#'   intercept random component and \eqn{\tau_1^2} is the variance of the treatment
#'   random effect.
#' @param rsq1 Estimate of variance explained by the level-1 (e.g., individual-level) covariates.
#' @param rsq2 Estimate of variance explained by the cluster-level covariates.
#' @param J Number of clusters. Determine \code{n} if \code{J} is specified.
#' @param n Cluster size. Determine \code{J} if \code{n} is specified.
#' @param K Number of cluster-level covariates.
#' @param P Proportion of the clusters that is treatment group.
#' @param alpha Significance level.
#' @return The expected confidence interval half width for two-level multisite randomized trials.
#' @export
#'

prec_msrt2 <- function(rho, omega, J, n, rsq1 = 0, rsq2 = 0,
                       K = 0, P = .5, alpha = .05) {
  df <- J - K - 1
  t_crit <- qt(1 - alpha / 2, df = df)
  t_crit * sqrt((rho * omega * (1 - rsq2) * P * (1 - P) * n +
                   (1 - rho) * (1 - rsq1)) /
                  (P * (1 - P) * J * n))
}
