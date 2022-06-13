#' Prior Distributions of Parameters
#'
#' \code{plot_prior()} generate plots of the prior distributions of effect size
#' and/or intraclass correlation.
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
#' @export

plot_prior <- function(d_est, d_sd, rho_est, rho_sd) {
  if (d_sd != 0) {
    pd <- ggplot2::ggplot(data.frame(d = c(d_est + 3 * d_sd, d_est - 3 * d_sd)),
                          ggplot2::aes(x = d)) +
      ggplot2::stat_function(fun = stats::dnorm, n = 101,
                             args = list(mean = d_est, sd = d_sd)) +
      ggplot2::labs(x = expression(delta), y = "Density")
  }
  if (rho_sd != 0) {
    shapes <- get_ab(rho_est, rho_sd)
    prho <- ggplot2::ggplot(data.frame(rho = c(rho_est + 4 * rho_sd, 0)),
                            ggplot2::aes(x = rho)) +
      ggplot2::stat_function(fun = stats::dbeta, n = 101,
                             args = list(shape1 = shapes[1], shape2 = shapes[2])) +
      ggplot2::labs(x = expression(rho), y = "Density")
  }
  if (d_sd != 0 & rho_sd == 0) {
    return(pd)
  } else if (d_sd == 0 & rho_sd != 0) {
    return(prho)
  } else if (d_sd != 0 & rho_sd != 0){
    return(list(delta = pd, rho = prho))
  } else {
    return("No prior distributions are specified. ")
  }
}
