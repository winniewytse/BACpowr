#' Determine Required Number of Clusters or Cluster Size for Two-Level CRTs
#'
#' \code{crtJn()} solves for the required number of clusters (J) or cluster size (n)
#' for a two-level CRT using the HCB approach. When certain uncertainty is specified,
#' this function determines the minimum required sample size that achieves the
#' desired level of expected power.
#'
#' @param d_est Effect size estimate.
#' @param d_sd Uncertainty in the effect size estimate.
#' @param rho_est Intraclass correlation estimate.
#' @param rho_sd Uncertainty in the intraclass correlation estimate.
#' @param r2_est Estimate of variance explained by the cluster-level covariates.
#' @param r2_sd Uncertainty in the variance explained by the cluster-level covariates.
#' @param J Specified number of clusters.
#' @param n Specified cluster size.
#' @param K Number of cluster-level covariates.
#' @param power Desired power level to achieve.
#' @param test One-tailed or two-tailed test.
#' @param plot Printing out a plot if it is TURE.
#' @param abs.tol Absolute tolerance. Defaults to \code{1e-10}.
#' @param x.tol X tolerance. Defaults to \code{1.5e-15}.
#' @param rel.tol Relative tolerance. Defaults to \code{1e-15}.
#' @param sing.tol Singular convergence tolerance. Defaults to \code{1e-20}.
#' @return The required J or n and a optionally plot that shows the power curve.
#' @export
#' @examples
#' crtJn(d_est = .5, d_sd = .2, rho_est = .1, rho_sd = .05, J = 30)
#' @seealso \url{https://winnie-wy-tse.shinyapps.io/hcb_shiny/}
crtJn <- function(d_est, d_sd, rho_est, rho_sd,
                  r2_est = 0, r2_sd = 0,
                  J = NULL, n = NULL, K = 0, power = .80,
                  test = "two-tailed", plot = FALSE,
                  abs.tol = 1e-10, x.tol = 1.5e-15,
                  rel.tol = 1e-15, sing.tol = 1e-20){
  ggplot2::theme_set(ggplot2::theme_bw())
  lossJ <- function(J) {
    sum((ep_crt(J = J, n = n, d_est = d_est, d_sd = d_sd,
                rho_est = rho_est, rho_sd = rho_sd,
                r2_est = r2_est, r2_sd = r2_sd,
                test = test) - power)^2)
  }
  lossn <- function(n) {
    sum((ep_crt(J = J, n = n, d_est = d_est, d_sd = d_sd,
                rho_est = rho_est, rho_sd = rho_sd,
                r2_est = r2_est, r2_sd = r2_sd,
                test = test) - power)^2)
  }
  if (!is.null(n)) {
    J <- stats::nlminb(start = 4, lossJ, lower = 1,
                       control = list(abs.tol = abs.tol, x.tol = x.tol,
                                      rel.tol = rel.tol, sing.tol = sing.tol))$par
  } else if (!is.null(J)) {
    n <- stats::nlminb(start = 0, lossn, lower = 1,
                       control = list(abs.tol = abs.tol, x.tol = x.tol,
                                      rel.tol = rel.tol, sing.tol = sing.tol))$par
  } else {
    n <- 1e10
    J <- stats::nlminb(start = c(4), lossJ, lower = c(1),
                       control = list(abs.tol = abs.tol, x.tol = x.tol,
                                      rel.tol = rel.tol, sing.tol = sing.tol))$par
    rm(n)
    n <- stats::nlminb(start = 0, lossn, lower = 1,
                       control = list(abs.tol = abs.tol, x.tol = x.tol,
                                      rel.tol = rel.tol, sing.tol = sing.tol))$par
  }
  if (plot) {
    p1 <- ggplot2::ggplot(data.frame(J = c(4, J + J/3)), ggplot2::aes(x = J)) +
      ggplot2::stat_function(fun = Vectorize(ep_crt, vectorize.args = "J"),
                             args = list(n = n, r2_est = r2_est, r2_sd = r2_sd,
                                         d_est = d_est, d_sd = d_sd,
                                         rho_est = rho_est, rho_sd = rho_sd),
                             n = 51) +
      ggplot2::geom_segment(x = J, xend = J, y = 0, yend = power,
                            linetype = "dashed", col = "red") +
      ggplot2::geom_segment(x = 0, xend = J, y = power, yend = power,
                            linetype = "dashed", col = "red") +
      ggplot2::labs(x = "Number of Clusters (J)", y = "Generalized Power")
    p2 <- ggplot2::ggplot(data.frame(n = c(1, n + n/3)), ggplot2::aes(x = n)) +
      ggplot2::stat_function(fun = Vectorize(ep_crt, vectorize.args = "n"),
                             args = list(J = J, r2_est = r2_est, r2_sd = r2_sd,
                                         d_est = d_est, d_sd = d_sd,
                                         rho_est = rho_est, rho_sd = rho_sd),
                             n = 51) +
      ggplot2::geom_segment(x = n, xend = n, y = 0, yend = power,
                            linetype = "dashed", col = "red") +
      ggplot2::geom_segment(x = 0, xend = n, y = power, yend = power,
                            linetype = "dashed", col = "red") +
      ggplot2::labs(x = "Cluster Size (n)", y = "Generalized Power")
    return(list(p1, p2, ceiling(cbind(J = J, n = n))))
  } else {
    return(ceiling(cbind(J = J, n = n)))
  }
}
