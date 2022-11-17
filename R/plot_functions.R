# This file contains plotting functions plot_prior(), plot_Jn() and plot_n().

#' Plot the Prior Distributions of Parameters
#'
#' \code{plot_prior()} generates plots of the prior distributions of effect size
#' and/or intraclass correlation.
#'
#' @param delta Effect size estimate, defined as
#'   \eqn{\delta = \frac{\gamma_{01}}{\tau^2 + \sigma^2}},
#'   where \eqn{\gamma_{01}} is the main effect of the treatment on the outcome,
#'   \eqn{\tau^2} is the variance of the cluster-specific random effect
#'   in the unconditional model (without covariates), and
#'   \eqn{\sigma^2} is the variance of the random error in the unconditional model.
#' @param delta_sd Uncertainty level of the effect size estimate.
#' @param rho Intraclass correlation estimate, defined as
#'   \eqn{\rho = \frac{\tau^2}{\tau^2 + \sigma^2}}, where \eqn{\tau^2} and
#'   \eqn{\sigma^2} are the variance components in the unconditional model.
#' @param rho_sd Uncertainty level of the intraclass correlation estimate.
#' @param omega Estimate of the treatment effect heterogeneity, defined as
#'   \eqn{\omega = \frac{\tau_1^2}{\tau_0^2}} where \eqn{\tau_0^2} is the
#'   variance of the intercept random component and \eqn{\tau_1^2} is the
#'   variance of the treatment random effect.
#' @param omega_sd Uncertainty level of the treatment effect heterogeneity estimate.
#' @export

plot_prior <- function(delta, delta_sd, rho, rho_sd, omega = NULL, omega_sd = NULL) {

  ggplot2::theme_set(ggplot2::theme_bw())

  if (is.null(omega_sd)) omega_sd <- 0
  if (is.null(rho_sd)) rho_sd <- 0
  if (sum(c(delta_sd, rho_sd, omega_sd) != 0) == 0) {
    stop("No prior distributions specified. ")
  }
  plots <- NULL
  if (delta_sd != 0) {
    pd <- ggplot2::ggplot(data.frame(d = c(delta + 3 * delta_sd, delta - 3 * delta_sd)),
                          ggplot2::aes(x = d)) +
      ggplot2::stat_function(fun = stats::dnorm, n = 101,
                             args = list(mean = delta, sd = delta_sd)) +
      ggplot2::labs(x = expression(delta), y = "Density")
    plots <- append(plots, list(delta = pd))
  }
  if (rho_sd != 0) {
    rho_ab <- beta_ab(rho, rho_sd)
    prho <- ggplot2::ggplot(data.frame(rho = c(rho + 4 * rho_sd, 0)),
                            ggplot2::aes(x = rho)) +
      ggplot2::stat_function(fun = stats::dbeta, n = 101,
                             args = list(shape1 = rho_ab[1], shape2 = rho_ab[2])) +
      ggplot2::labs(x = expression(rho), y = "Density")
    plots <- append(plots, list(rho = prho))
  }
  if (omega_sd != 0) {
    omega_ab <- gamma_ab(omega, omega_sd)
    pomega <- ggplot2::ggplot(data.frame(omega = c(omega + 4 * omega_sd, 0)),
                              ggplot2::aes(x = omega)) +
      ggplot2::stat_function(fun = stats::dgamma, n = 101,
                             args = list(shape = omega_ab[1], rate = omega_ab[2])) +
      ggplot2::labs(x = expression(omega), y = "Density")
    plots <- append(plots, list(omega = pomega))
  }
  return(plots)
}



# Plot n

plot_n <- function(n, delta, delta_sd, power = .8, alpha = .05,
                   test = "two.sided", ep = NULL, al = NULL, smooth = 51) {

  if (is.null(al)) {
    fn <- ep_indp_t
    if (delta_sd == 0) {
      ylab <- "Classical Power"
    } else {
      ylab <- "Expected Power"
    }
    yval <- power
  } else {
    fn <- al_indp_t
    ylab <- "Assurance Level"
    yval <- al
  }

  ggplot2::theme_set(ggplot2::theme_bw())

  pn <- ggplot2::ggplot(data.frame(n = c(3, n + n / 5)), ggplot2::aes(x = n)) +
    ggplot2::stat_function(fun = Vectorize(fn, vectorize.args = "n1"),
                           args = list(delta = delta, delta_sd = delta_sd,
                                       n2 = NULL, alpha = alpha, power = power,
                                       test = test),
                           n = smooth) +
    ggplot2::geom_segment(x = n, xend = n, y = 0, yend = yval,
                          linetype = "dashed", col = "red") +
    ggplot2::geom_segment(x = 0, xend = n, y = yval, yend = yval,
                          linetype = "dashed", col = "red") +
    ggplot2::labs(x = "Group Size (n)", y = ylab)

  return(list(n = pn))
}


# Plot Jn

plot_Jn <- function(J, n, delta, delta_sd, rho, rho_sd,
                    omega = NULL, omega_sd =  NULL, rsq1 = 0, rsq2 = 0,
                    K = 0, P = .5, power = .8, alpha = .05,
                    test = "two.sided", ep = NULL, al = NULL, smooth = 51) {

  if (is.null(omega_sd)) {
    omega_sd <- 0
    args <- list(delta = delta, delta_sd = delta_sd, rho = rho, rho_sd = rho_sd,
                 rsq2 = rsq2, K = K, P = P, power = power, alpha = alpha, test = test)
    if (is.null(al)) {
      fn <- ep_crt2
    } else {
      fn <- al_crt2
    }
  } else {
    args <- list(delta = delta, delta_sd = delta_sd, rho = rho, rho_sd = rho_sd,
                 omega = omega, omega_sd = omega_sd, rsq1 = rsq1, rsq2 = rsq2, K = K,
                 P = P, power = power, alpha = alpha, test = test)
    if (is.null(al)) {
      fn <- ep_msrt2
    } else {
      fn <- al_msrt2
    }
  }

  if (is.null(al)) {
    if (delta_sd == 0 & rho_sd == 0 & omega_sd == 0) {
      y_lab <- "Classical Power"
    } else {
      y_lab <- "Expected Power"
    }
    yval <- power
  } else {
    y_lab <- "Assurance Level"
    yval <- al
  }

  ggplot2::theme_set(ggplot2::theme_bw())

  pJ <- ggplot2::ggplot(data.frame(J = c(J - J / 2, J + J / 5)),
                        ggplot2::aes(x = J)) +
    ggplot2::stat_function(fun = Vectorize(fn, vectorize.args = "J"),
                           args = append(args, list(n = n)), n = smooth) +
    ggplot2::geom_segment(x = J, xend = J, y = 0, yend = yval,
                          linetype = "dashed", col = "red") +
    ggplot2::geom_segment(x = 0, xend = J, y = yval, yend = yval,
                          linetype = "dashed", col = "red") +
    ggplot2::labs(x = "Number of Clusters (J)", y = y_lab)

  pn <- ggplot2::ggplot(data.frame(n = c(1, n + n / 5)), ggplot2::aes(x = n)) +
    ggplot2::stat_function(fun = Vectorize(fn, vectorize.args = "n"),
                           args = append(args, list(J = J)), n = smooth) +
    ggplot2::geom_segment(x = n, xend = n, y = 0, yend = yval,
                          linetype = "dashed", col = "red") +
    ggplot2::geom_segment(x = 0, xend = n, y = yval, yend = yval,
                          linetype = "dashed", col = "red") +
    ggplot2::labs(x = "Cluster Size (n)", y = y_lab)

  return(list(J = pJ, n = pn))
}
