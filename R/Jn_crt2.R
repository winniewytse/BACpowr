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
#' @return The required J or n and a optionally plot that shows the power curve.
#' @export
#' @examples
#' Jn_crt2(d_est = .5, d_sd = .2, rho_est = .1, rho_sd = .05, J = 30)
#' @seealso \url{https://winnie-wy-tse.shinyapps.io/hcb_shiny/}

Jn_con_crt2 <- function(d_est, rho_est, r2_est = 0, J = NULL, n = NULL,
                        K = 0, power = .80, test = "two-tailed") {
  lossJ <- function(J) {
    sum((pow_crt2(J = J, n = n, d_est = d_est, rho_est = rho_est,
                  r2_est = r2_est, test = test) - power)^2)
  }
  lossn <- function(n) {
    sum((pow_crt2(J = J, n = n, d_est = d_est, rho_est = rho_est,
                  r2_est = r2_est, test = test) - power)^2)
  }
  if (is.null(J)) {
    lbfgsb <- optim(K + 2 + 1, lossJ, lower = K + 2 + 1, upper = Inf,
                    method = "L-BFGS-B")
  } else if (is.null(n)) {
    lbfgsb <- optim(1, lossn, lower = 1, upper = Inf, method = "L-BFGS-B")
  }
  lbfgsb$par
}

Jn_crt2 <- function(d_est, d_sd, rho_est, rho_sd, r2_est = 0, r2_sd = 0,
                    J = NULL, n = NULL, K = 0, power = .80, al = NULL,
                    maxiter = 100,
                    test = "two-tailed", plot = FALSE){

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
    if (is.null(J)) {
      minJ <- K + 2 + 1
    }
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
    if (is.null(J)) {
      # solve minimum J for a nonzero assurance level
      # to avoid being stuck at local minimum
      minJ <- Jn_con_crt2(d_est, rho_est, r2_est, J, n, K, power, test)[1]
    } else {
      minJ <- NULL
    }
  }

  if (is.null(J)) {

    brent <- optim(minJ, lossJ, lower = minJ, upper = 1e6, method = "Brent")
    # if L-BFGS-B does not converge, try using PORT routines
    if (brent$value > 1e-3) {
      port <- nlminb(minJ, lossJ, lower = minJ)
      # J <- optimize(lossJ, c(minJ, 1e6))$minimum
      if (port$objective > 1e-3) {
        lbfgsb <- optim(minJ, lossJ, lower = minJ, upper = Inf, method = "L-BFGS-B")
        if (lbfgsb$value > 1e-3) {
          J <- lbfgsb$par
          warning(paste0("The algorithm fails to converge for the specified priors. ",
                         "There may not exist a solution for the desired expected ",
                         "power or assurance level. ",
                         "Please consider some lower power/assurance level. "))
        } else {
          J <- lbfgsb$par
        }
      } else {
        J <- port$par
      }
    } else {
      J <- brent$par
    }

  } else if (is.null(n)) {

    lbfgsb <- optim(1, lossn, lower = 1, upper = Inf, method = "L-BFGS-B")
    # if L-BFGS-B does not converge, try using PORT routines
    if (lbfgsb$value > 1e-3) {
      port <- nlminb(1, lossn, lower = 1)
      # if nlminb fails as well
      if (port$object > 1e-3) {
        stop(paste0("The algorithm fails to converge due to too few J ",
                    "for the specified priors. \n",
                    "Please consider raising J."))
      } else {
        n <- port$par
      }
    } else {
      n <- lbfgsb$par
    }
  } else {
    n <- 1e10
    J <- stats::nlminb(start = minJ, lossJ, lower = 1)$par
    rm(n)
    n <- stats::nlminb(start = 0, lossn, lower = 1)$par
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
