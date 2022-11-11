#' Classical Statistical Power for Two-Level CRTs
#'
#' \code{pow_crt2()} computes classical statistical power for a two-level CRT.
#'
#' @param delta Effect size estimate, defined as
#'   \eqn{\delta = \frac{\gamma_{01}}{\tau^2 + \sigma^2}},
#'   where \eqn{\gamma_{01}} is the main effect of the treatment on the outcome,
#'   \eqn{\tau^2} is the variance of the cluster-specific random effect
#'   in the unconditional model (without covariates), and
#'   \eqn{\sigma^2} is the variance of the random error in the unconditional model.
#' @param rho Intraclass correlation estimate, defined as
#'   \eqn{\rho = \frac{\tau^2}{\tau^2 + \sigma^2}}, where \eqn{\tau^2} and
#'   \eqn{\sigma^2} are the variance components in the unconditional model.
#' @param rsq2 Estimate of variance explained by the cluster-level covariates.
#' @param J Number of clusters.
#' @param n Cluster size.
#' @param K Number of cluster-level covariates.
#' @param P Proportion of the clusters that is treatment group.
#' @param alpha Type I error rate. Default to be \code{.05}.
#' @param test One-sided or two-sided test. Options are either "one.sided" or
#' "two.sided".
#' @return Classical statistical power for a two-level CRT design with J
#'  clusters each has n observations.
#' @export
#' @examples
#' pow_crt2(J = 30, n = 100, delta = .5, rho = .1)
#' pow_crt2(J = 30, n = 100, delta = .5, rho = .1, rsq2 = .3)
pow_crt2 <- function(J, n, delta, rho, rsq2 = 0,
                     K = 0, P = .5, alpha = .05,
                     test = "two.sided") {

  if (J <= K + 2) stop(paste0("J needs to be larger than the number of parameters. ",
                              "Please increase J."))
  df <- J - K - 2
  ncp <- delta * sqrt(J * n * P * (1 - P) / (1 + (n * (1 - rsq2) - 1) * rho))

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

