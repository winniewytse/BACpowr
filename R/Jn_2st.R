#' Determine Number of Clusters or Cluster Size for Two-Level Multisite Randomzied Trials
#'
#' \code{Jn_2st()} solves for the required group size for a two-sample t-test.
#' When the uncertainty level of the effect size is specified, this function
#' determines the sample size requisite that achieves the desired expected power
#' or assurance level. Otherwise, this function determines the sample size requisite
#' that achieves the desired classical power.
#'
#' @param d_est Effect size estimate, defined as
#'   \eqn{d = \frac{\bar x_2 - \bar x_1}
#'   {\sqrt{\frac{(n_1 - 1)s_1^2 + (n_2 - 1)s_2^2}{n_1 + n_2 - 2}}}}
#' @param d_sd Uncertainty level of the effect size estimate.
#' @param power Desired level of statistical power.
#' @param alpha Type I error rate. Default to be \code{.05}.
#' @param ep Desired expected power to achieve. If neither \code{ep} nor
#'   \code{al} is specified, \code{ep} = \code{power}. An 80% expected power
#'   indicates that the mean or average power value is 80% over
#'   the specified uncertainty.
#' @param al Assurance level to achieve. An 80% assurance level indicates 80% of the
#'   power values are above the desired statistical power over the specified uncertainty.
#'   Default to be \code{.6}.
#' @param test One-sided or two-sided test. Options are either "one.sided" or "two.sided".
#' @param plot If TRUE, plots of J and n against the expected power or assurance level
#'   will be printed.
#' @return The required group size for a two-sample t-test.
#' @import stats
#' @export
#' @examples
#' Jn_2st(d_est = .5, d_sd = .1)
#' @seealso \url{https://winnie-wy-tse.shinyapps.io/hcb_shiny/}

Jn_2st <- function(d_est, d_sd,
                   alpha = .05, power = .8, ep = NULL, al = NULL,
                   test = "two.sided", plot = FALSE) {

  ggplot2::theme_set(ggplot2::theme_bw())

  if (is.null(ep) & is.null(al)) ep <- power

  if (is.null(al)) { # solve with the expected power
    criteria <- ep_2st
    target <- ep
  } else { # solve with the assurance level
    criteria <- al_2st
    target <- al
  }

  loss <- function(n) {
    criteria(d_est = d_est, d_sd = d_sd, n1 = n, n2 = n,
             alpha = alpha, power = power, test = test) - target
  }
  n <- try(stats::uniroot(loss, c(3, 1e8))$root, silent = TRUE)
  # if root-finding method fails, try optimization methods
  if (class(n) == "try-error") {
    loss <- function(n) {
      (criteria(d_est = d_est, d_sd = d_sd, n1 = n, n2 = n,
                alpha = alpha, test = test) - target)^2
    }
    n <- Jn_optimize(start = 3, loss = loss, lower = 1, upper = 1e6)
  }


  if (plot) {
    n_plot <- plot_n(n = n, d_est = d_est, d_sd = d_sd,
                     power = power, alpha = alpha,
                     ep = ep, al = al)
    prior_plot <- plot_prior(d_est = d_est, d_sd = d_sd,
                             rho_est = NULL, rho_sd = NULL)

    return(list(n_plot = n_plot, prior_plot = prior_plot,
                n = ceiling(n)))

  } else {
    return(ceiling(n))
  }
}
