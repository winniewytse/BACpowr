#' Determine Cluster Size (n) for Independent Sample T-Tests
#'
#' \code{n_indp_t()} solves for the required cluster size for a two-sample t-test.
#' When the uncertainty level of the effect size is specified, determines the
#' sample size requisite that achieves the desired expected power or assurance
#' level. Otherwise, determines the sample size requisite that achieves the
#' desired classical power.
#'
#' @param delta Effect size estimate, defined as
#'   \eqn{d = \frac{\bar x_2 - \bar x_1}
#'   {\sqrt{\frac{(n_1 - 1)s_1^2 + (n_2 - 1)s_2^2}{n_1 + n_2 - 2}}}}
#' @param delta_sd Uncertainty level of the effect size estimate.
#' @param power Desired statistical power level.
#' @param alpha Type I error rate. Defaults to \code{.05}.
#' @param ep Desired expected power (EP). For example, an 80% EP
#'  indicates that the mean or average power value is 80% over the specified
#'  uncertainty. If neither \code{ep} nor \code{al} is specified,
#'  \code{ep} = \code{power}.
#' @param al Desired assurance level (AL). For example, an 80% AL indicates
#'  that 80% of the power values are above the desired statistical power over
#'  the specified uncertainty.
#' @param test Whether a one-sided or two-sided test should be performed.
#'   Defaults to "two-sided".
#' @param plot If TRUE, plots of J and n against EP or AL will be printed.
#' @return Required group size for a two-sample t-test.
#' @import stats
#' @export
#' @examples
#' n_indp_t(delta = .5, delta_sd = .1)
#' @seealso \url{https://winnie-wy-tse.shinyapps.io/hcb_shiny/}

n_indp_t <- function(delta, delta_sd, alpha = .05, power = .8, ep = NULL,
                  al = NULL, test = "two.sided", plot = FALSE) {

  if (is.null(ep) & is.null(al)) ep <- power

  if (is.null(al)) { # solve with the expected power
    criteria <- ep_indp_t; target <- ep
  } else { # solve with the assurance level
    criteria <- al_indp_t; target <- al
  }

 # params <- list(delta = delta, delta_sd = delta_sd, n1 = n, n2 = n,
  #               alpha = alpha, test = test)

  loss <- function(n) {
   # do.call(criteria, append(list(power=power), params)) - target
    criteria(delta = delta, delta_sd = delta_sd, n1 = n, n2 = n,
             alpha = alpha, power = power, test = test) - target
    }
  n <- try(stats::uniroot(loss, c(3, 1e8))$root, silent = TRUE)

  # if root-finding method fails, try optimization methods
  if (class(n) == "try-error") {
    loss <- function(n) {
   #   (do.call(criteria, append(list(power=power), params)) - target)^2
      (criteria(delta = delta, delta_sd = delta_sd, n1 = n, n2 = n,
                alpha = alpha, test = test) - target)^2
    }
    n <- optimize_Jn(start = 3, loss = loss, lower = 1, upper = 1e6,
                     solve = "n")
  }

  if (plot) {
    n_plot <- plot_n(n = n, delta = delta, delta_sd = delta_sd, power = power,
                     alpha = alpha, ep = ep, al = al)
    if (delta_sd == 0) {
      prior_plot <- NULL
    } else {
      prior_plot <- plot_prior(delta = delta, delta_sd = delta_sd,
                               rho = NULL, rho_sd = NULL)
    }
    return(list(n_plot = n_plot, prior_plot = prior_plot, n = ceiling(n)))

  } else {
    return(ceiling(n))
  }
}
