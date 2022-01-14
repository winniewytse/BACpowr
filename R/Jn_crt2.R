#' Determine Required Number of Clusters or Cluster Size for Two-Level CRTs
#'
#' \code{Jn_crt2()} solves for the required number of clusters (J) or cluster size (n)
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
#' @param criteria Solve the required J or n based on the desired assurance level
#' or expected power
#' @param test One-tailed or two-tailed test.
#' @param plot Printing out a plot if it is TURE.
#' @param abs.tol Absolute tolerance. Defaults to \code{1e-10}.
#' @param x.tol X tolerance. Defaults to \code{1.5e-15}.
#' @param rel.tol Relative tolerance. Defaults to \code{1e-15}.
#' @param sing.tol Singular convergence tolerance. Defaults to \code{1e-20}.
#' @return The required J or n and a optionally plot that shows the power curve.
#' @export
#' @examples
#' Jn_crt2(d_est = .5, d_sd = .2, rho_est = .1, rho_sd = .05, J = 30)
#' @seealso \url{https://winnie-wy-tse.shinyapps.io/hcb_shiny/}
Jn_crt2 <- function(d_est, d_sd, rho_est, rho_sd,
                    r2_est = 0, r2_sd = 0,
                    J = NULL, n = NULL, K = 0,
                    power = .80, al = NULL,
                    test = "two-tailed", plot = FALSE,
                    abs.tol = 1e-10, x.tol = 1.5e-15,
                    rel.tol = 1e-15, sing.tol = 1e-20){

  if (is.null(al)) {
    lossJ <- function(J) {
      sum((ep_crt2(J = J, n = n, d_est = d_est, d_sd = d_sd,
                   rho_est = rho_est, rho_sd = rho_sd,
                   r2_est = r2_est, r2_sd = r2_sd,
                   test = test) - power)^2)
    }
    lossn <- function(n) {
      sum((ep_crt2(J = J, n = n, d_est = d_est, d_sd = d_sd,
                   rho_est = rho_est, rho_sd = rho_sd,
                   r2_est = r2_est, r2_sd = r2_sd,
                   test = test) - power)^2)
    }
    minJ <- 4
  } else if (!is.null(al)) {
    lossJ <- function(J) {
      sum(al_crt2(J = J, n = n, d_est = d_est, d_sd = d_sd,
                  rho_est = rho_est, rho_sd = rho_sd,
                  r2_est = r2_est, r2_sd = r2_sd,
                  test = test) - al)^2
    }
    lossn <- function(n) {
      sum((al_crt2(J = J, n = n, d_est = d_est, d_sd = d_sd,
                   rho_est = rho_est, rho_sd = rho_sd,
                   r2_est = r2_est, r2_sd = r2_sd,
                   test = test) - al)^2)
    }
    # solve minimum J for a nonzero assurance level
    # to avoid being stuck at local minimum
    minJ <- 4
    a <- al_crt2(J = minJ, n = n, d_est = d_est, d_sd = d_sd,
                 rho_est = rho_est, rho_sd = rho_sd,
                 r2_est = r2_est, r2_sd = r2_sd,
                 test = test)
    while (a < 1e-4) {
      minJ <- minJ + 1
      a <- al_crt2(J = minJ, n = n, d_est = d_est, d_sd = d_sd,
                   rho_est = rho_est, rho_sd = rho_sd,
                   r2_est = r2_est, r2_sd = r2_sd,
                   test = test)
    }
  }

  if (!is.null(n)) {
    output <- optim(minJ, lossJ, lower = minJ, upper = Inf, method = "L-BFGS-B")
    if (output$value > 1e-3) {
      J <- optimize(lossJ, c(minJ, 1e6))$minimum
    } else {
      J <- output$par
    }
  } else if (!is.null(J)) {
    n <- optim(1, lossn, lower = 1, upper = Inf, method = "L-BFGS-B")$par
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

  if (J >= 9e5) warning(paste0("The minimum J requisite may be unreasonably large. ",
                               "Please check if the priors are correctly specified."))
  if (n > 9e5) warning(paste0("The minimum n requisite may be unreasonably large. ",
                              "Please consider raising J. "))

  if (plot) {
    ggplot2::theme_set(ggplot2::theme_bw())

    if (is.null(al)) {
      plots <- Jn_plot(J = J, n = n, d_est = d_est, d_sd = d_sd,
                       rho_est = rho_est, rho_sd = rho_sd,
                       r2_est = 0, r2_sd = 0, power = power, al = NULL)
    } else {
      plots <- Jn_plot(J = J, n = n, d_est = d_est, d_sd = d_sd,
                       rho_est = rho_est, rho_sd = rho_sd,
                       r2_est = 0, r2_sd = 0, power = power,
                       al = al, minJ = minJ)
    }

    if (J >= 9e5)
      warning(paste0("Plots may be unreliable."))

    return(list(plots, ceiling(cbind(J = J, n = n))))

  } else {
    return(ceiling(cbind(J = J, n = n)))
  }
}
