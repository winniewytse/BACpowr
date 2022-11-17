#' Determine Number of Clusters or Cluster Size for Two-Level CRTs
#'
#' \code{Jn_crt2()} uses the HCB approach to solve for the minimum required
#' number of clusters (J) or cluster size (n) that would achieve a desired
#' expected power or assurance level for a two-level CRT.
#'
#' @param delta Effect size estimate, defined as
#'   \eqn{\delta = \frac{\gamma_{01}}{\tau^2 + \sigma^2}}, where
#'   \eqn{\gamma_{01}} is the main effect of the treatment on the outcome,
#'   \eqn{\tau^2} is the variance of the cluster-specific random effect in the
#'   unconditional model (without covariates), and \eqn{\sigma^2} is the
#'   variance of the random error in the unconditional model.
#' @param delta_sd Uncertainty level of the effect size estimate.
#' @param rho Intraclass correlation (ICC) estimate, defined as
#'   \eqn{\rho = \frac{\tau^2}{\tau^2 + \sigma^2}}, where \eqn{\tau^2} and
#'   \eqn{\sigma^2} are the variance components in the unconditional model.
#' @param rho_sd Uncertainty level of the ICC estimate.
#' @param rsq2 Estimate of variance explained by the cluster-level covariates.
#'   \code{0} by default.
#' @param J Number of clusters. If \code{J} is specified, \code{Jn_crt2()}
#'   determines \code{n} for the specified expected power/assurance level.
#' @param n Cluster size. If \code{n} is specified, \code{Jn_crt2()}
#'   determines \code{J} for the specified expected power/assurance level.
#' @param K Number of cluster-level covariates. \code{0} by default.
#' @param P Proportion of the clusters that are treatment groups. \code{.5} by
#'   default.
#' @param power Desired level of statistical power. \code{.8} by default.
#' @param alpha Type I error rate. \code{.05} by default.
#' @param ep Desired expected power (EP). If neither \code{ep} nor \code{al} is
#'   specified, \code{ep} = \code{power}. For example, an 80% EP indicates that
#'   the mean or average power value is 80% over the specified uncertainty.
#' @param al Desired assurance level (AL). For example, an 80% AL indicates 80%
#'   of the power values are above the desired statistical power over the
#'   specified uncertainty.
#' @param test Whether a one-sided or two-sided test should be performed.
#' @param plot (optional) whether plots of J and n against expected power or
#'   assurance level should be returned.
#' @return A 1 x 2 array containing the J and n values that together achieve the
#'   desired expected power or assurance level. If \code{plot = TRUE}, also
#'   returns plots.
#' @import stats
#' @export
#' @examples
#' Jn_crt2(delta = .5, delta_sd = .2, rho = .1, rho_sd = .05, J = 30)
#' @seealso \url{https://winnie-wy-tse.shinyapps.io/hcb_shiny/}

Jn_crt2 <- function(delta, delta_sd, rho, rho_sd, rsq2 = 0, J = NULL,
                    n = NULL, K = 0, P = .5, alpha = .05, power = .8,
                    ep = NULL, al = NULL, test = "two.sided",
                    plot = FALSE) {

  if (is.null(J) & is.null(n)) stop(paste0("Please specify either n or J."))

  # If neither EP nor AL is specified, set EP equal to power to solve for power.
  if (is.null(ep) & is.null(al)) {ep <- power}

  # If both EP and AL were specified, solve sample size for desired AL.
  if (!is.null(ep) & !is.null(al)) {
    ep <- NULL
  }

  # As a starting point, compute J and n using the conventional approach.
  Jn_conv <- Jn_crt2_c(delta = delta, rho = rho, rsq2 = rsq2,
                       J = J, n = n, K = K, P = P, alpha = alpha, power = power,
                       test = test)

  # If uncertainty is set to 0 for effect size and ICC estimates, return J and n
  # values computed using the conventional approach and plots if plot == TRUE.
  if (delta_sd == 0 & rho_sd == 0) {
    if (plot) {
      Jn_plots <- plot_Jn(J = Jn_conv[1], n = Jn_conv[2], delta = delta,
                          delta_sd = delta_sd, rho = rho, rho_sd = rho_sd,
                          rsq2 = rsq2, K = K, P = P, power = power,
                          alpha = alpha, ep = ep, al = al)
      return(list(Jn_plots = Jn_plots, Jn = ceiling(Jn_conv)))
    } else {
      return(ceiling(Jn_conv))
    }
  }

  # If AL is not specified, solve with EP (target) using ep_crt2().
  if (is.null(al) & !is.null(ep)) {
    criteria <- ep_crt2; target <- ep
  }
  # If EP is not specified, solve with AL (target) using al_crt2().
  if (is.null(ep) & !is.null(al)) {
    criteria <- al_crt2; target <- al
  }

  params <- list(delta = delta, delta_sd = delta_sd, rho = rho, rho_sd = rho_sd,
                 rsq2 = rsq2, K = K, P = P, power = power, alpha = alpha)

  # Define a loss function for J or n, attempt to optimize using uniroot in the
  # specified internal, and try other optimization methods if root-finding fails.
  if (is.null(J)) { # solve for J
    loss <- function(J) {
      do.call(criteria, append(list(J = J, n = n), params)) - target
    }
    min_j <- if (is.null(al) & !is.null(ep)) (K + 2 + 1) else Jn_conv[1]

    J <- try(stats::uniroot(loss, interval = c(min_j, 1e8))$root, silent = TRUE)

    # if root-finding method fails, try optimization methods
    if (class(J) == "try-error") {
      loss <- function(J) {
        (do.call(criteria, append(list(J = J, n = n), params)) - target)^2
      }
      J <- optimize_Jn(start = min_j, loss = loss, lower = K + 3, upper = 1e6,
                       solve = "J")
    }
  } else { # solve n
    loss <- function(n) {
      do.call(criteria, append(list(J = J, n = n), params)) - target
    }
    min_n <- 1

    n <- try(stats::uniroot(loss, interval = c(min_n, 1e8))$root, silent = TRUE)
    if (class(n) == "try-error") {
      loss <- function(n) {
        (do.call(criteria, append(list(J = J, n = n), params)) - target)^2
      }
      n <- optimize_Jn(start = min_n, loss = loss, lower = 1, upper = Inf,
                       solve = "n")
    }
  }
  if (J >= 9e5) warning(paste0("Results may be unreliable due to convergence
                                 issues. Try increasing the cluster size (n),
                                 using smaller uncertainty (delta_sd or rho_sd),
                                 or decreasing EP or AL."))

  if (plot) {
    prior_plots <- plot_prior(delta = delta, delta_sd = delta_sd, rho = rho,
                              rho_sd = rho_sd)
    Jn_plots <- plot_Jn(J = J, n = n, delta = delta, delta_sd = delta_sd,
                        rho = rho, rho_sd = rho_sd, rsq2 = rsq2, K = K,
                        P = P, power = power, alpha = alpha, ep = ep, al = al)

    return(list(Jn_plots = Jn_plots, prior_plots = prior_plots,
                Jn = ceiling(cbind(J = J, n = n))))
  } else {
    return(ceiling(cbind(J = J, n = n)))
  }
}
