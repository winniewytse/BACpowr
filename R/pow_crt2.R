#' Statistical Power for Two-Level CRTs
#'
#' \code{pow_crt2()} computes the statistical power for a two-level CRT.
#'
#' @param d_est Effect size estimate.
#' @param rho_est Intraclass correlation estimate.
#' @param r2_est Estimate of variance explained by the cluster-level covariates.
#' @param J Specified number of clusters.
#' @param n Specified cluster size.
#' @param K Number of cluster-level covariates.
#' @param test One-tailed or two-tailed test.
#' @return Statistical power given certain J and n for cluster-randomized trials.
#' @export
#' @examples
#' pow_crt2(J = 30, n = 100, d_est = .5, rho_est = .1)
#' pow_crt2(J = 30, n = 100, d_est = .5, rho_est = .1, r2_est = .3)
pow_crt2 <- function(J, n, d_est, rho_est, r2_est = 0, K = 0, P = .5,
                     test = "two-tailed") {
  if (J <= K + 2) stop(paste0("J needs to be larger than the number of parameters. ",
                              "Please increase J."))
  df <- J - K - 2
  if (test == "two-tailed") {
    cv <- stats::qt(.975, df)
    ncp <- d_est * sqrt(J * n * P * (1 - P) / (1 + (n * (1 - r2_est) - 1) * rho_est))
    pow <- stats::pt(cv, df = df, ncp = ncp, lower.tail = FALSE) +
      stats::pt(-cv, df = df, ncp = ncp, lower.tail = TRUE)
  } else if (test == "one-tailed") {
    cv <- stats::qt(.95, df)
    ncp <- d_est * sqrt(J * n * P * (1 - P) / (1 + (n * (1 - r2_est) - 1) * rho_est))
    pow <- stats::pt(cv, df = df, ncp = ncp, lower.tail = FALSE)
  }
  return(pow)
}

