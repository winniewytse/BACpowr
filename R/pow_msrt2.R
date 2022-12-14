#' Classical Statistical Power for Two-Level MSRTs
#'
#' \code{pow_msrt2()} computes the classical statistical power for a two-level
#'  multisite randomized trial (MSRT).
#'
#' @param delta Effect size estimate, defined as
#'  \eqn{\delta = \frac{\gamma_{01}}{\tau^2 + \sigma^2}}, where
#'  \eqn{\gamma_{01}} is the main effect of the treatment on the outcome,
#'  \eqn{\tau^2} is the variance of the cluster-specific random effect in the
#'  unconditional model (without covariates), and \eqn{\sigma^2} is the
#'  variance of the random error in the unconditional model.
#' @param rho Intraclass correlation (ICC) estimate, defined as
#'  \eqn{\rho = \frac{\tau^2}{\tau^2 + \sigma^2}}, where \eqn{\tau^2} and
#'  \eqn{\sigma^2} are the variance components in the unconditional model.
#' @param omega Estimate of the treatment effect heterogeneity, defined as
#'  \eqn{\omega = \frac{\tau_1^2}{\tau_0^2}} where \eqn{\tau_0^2} is the
#'  variance of the intercept random component and \eqn{\tau_1^2} is the
#'  variance of the treatment random effect.
#' @param rsq1 Estimate of variance explained by the level-1 (individual-level)
#'  covariates. Defaults to \code{0}.
#' @param rsq2 Estimate of variance explained by the level-2 (cluster-level)
#'  covariates. Defaults to \code{0}.
#' @param J Number of clusters.
#' @param n Cluster size.
#' @param K Number of cluster-level covariates. Defaults to \code{0}.
#' @param P Proportion of the sample that is treatment group. \code{.5} by
#'   default.
#' @param alpha Type I error rate. Defaults to \code{.05}.
#' @param test Whether a one-sided or two-sided test should be performed.
#'   Defaults to "two-sided".
#' @return Classical statistical power for a two-level MSRT design.
#' @export
#' @examples
#' pow_msrt2(J = 30, n = 15, delta = .3, rho = .2, omega = .3)
#' pow_msrt2(J = 10, n = 80, delta = .3, rho = .2, omega = .2, P = .3)
pow_msrt2 <- function(J, n, delta, rho, omega, rsq1 = 0, rsq2 = 0, K = 0,
                      P = .5, alpha = .05, test = "two.sided") {
  if (J <= K + 1) stop(paste0("J needs to be larger than the number of parameters. ",
                              "Please increase J."))
  # if (rho == 1 & omega == 0)
  #   stop("The denominator of the ncp is 0 when rho = 1 and omega = 0.")

  df <- J - K - 1
  # non-centrality parameter (ncp)
  ncp <- delta * sqrt(P * (1 - P) * J * n /
                        (rho * omega * (1 - rsq2) * P * (1 - P) * n +
                           (1 - rho) * (1 - rsq1)))
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
