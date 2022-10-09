#' Determine Number of Clusters or Cluster Size for Two-Level CRTs
#'
#' \code{Jn_crt2()} uses the HCB approach to solve for the minimum required
#' number of clusters (J) or cluster size (n) that would achieve a desired
#' expected power or assurance level for a two-level CRT.
#'
#' @param d_est Effect size estimate, defined as
#'   \eqn{\delta = \frac{\gamma_{01}}{\tau^2 + \sigma^2}}, where
#'   \eqn{\gamma_{01}} is the main effect of the treatment on the outcome,
#'   \eqn{\tau^2} is the variance of the cluster-specific random effect in the
#'   unconditional model (without covariates), and \eqn{\sigma^2} is the
#'   variance of the random error in the unconditional model.
#' @param d_sd Uncertainty level of the effect size estimate.
#' @param rho_est Intraclass correlation (ICC) estimate, defined as
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
#' Jn_crt2(d_est = .5, d_sd = .2, rho_est = .1, rho_sd = .05, J = 30)
#' @seealso \url{https://winnie-wy-tse.shinyapps.io/hcb_shiny/}

Jn_crt2 <- function(d_est, d_sd, rho_est, rho_sd, rsq2 = 0,
                    J = NULL, n = NULL, K = 0, P = .5,
                    alpha = .05, power = .8, ep = NULL, al = NULL,
                    test = "two.sided", reparameterize = FALSE,
                    plot = FALSE) {

  if (is.null(ep) & is.null(al)) ep <- power

  # Compute J / n using the conventional approach as a starting point.
  Jn_crt <- Jn_crt2_c(d_est = d_est, rho_est = rho_est, rsq2 = rsq2,
                      J = J, n = n, K = K, P = P,
                      alpha = alpha, power = power, test = test,
                      reparameterize = reparameterize)

  # If the uncertainty for the effect size estimate and the ICC estimate are both
  # 0, return the J / n combination computed using the conventional approach
  # and plots if plot = TRUE.
  if (d_sd == 0 & rho_sd == 0) {
    if (plot) {
      Jn_plots <- plot_Jn(J = Jn_crt[1], n = Jn_crt[2],
                          d_est = d_est, d_sd = d_sd,
                          rho_est = rho_est, rho_sd = rho_sd,
                          rsq2 = rsq2, K = K, P = P,
                          power = power, alpha = alpha, ep = ep, al = al)
      return(list(Jn_plots = Jn_plots, Jn = ceiling(Jn_crt)))
    } else {
      return(ceiling(Jn_crt))
    }
  }

  # If assurance level is not specified, Jn_crt2() will solve with expected
  # power (target) using function ep_crt2() (criteria).
  if (is.null(al)) {
    criteria <- ep_crt2
    target <- ep
    if (is.null(J)) {
      min <- K + 2 + 1 # degrees of freedom (DF = J-K-2)
    }
  } else {
  # If expected power is not specified, Jn_crt2() will solve with assurance
  # level (target) using function al_crt2() (criteria).
    criteria <- al_crt2
    target <- al
    if (is.null(J)) {
      # set a higher min J to avoid being stuck at the local minimum
      min <- Jn_crt[1] # Based on the conventional approach
    }
  }


  if (is.null(J)) { # solve J
    loss <- function(J) {
      criteria(J = J, n = n, d_est = d_est, d_sd = d_sd, rho_est = rho_est,
               rho_sd = rho_sd, rsq2 = rsq2, K = K, P = P, power = power,
               alpha = alpha, test = test, reparameterize = reparameterize) - target
    }
    J <- try(stats::uniroot(loss, c(min, 1e8))$root, silent = TRUE)
    # if root-finding method fails, try optimization methods
    if (class(J) == "try-error") {
      loss <- function(J) {
        (criteria(J = J, n = n, d_est = d_est, d_sd = d_sd, rho_est = rho_est,
                  rho_sd = rho_sd, rsq2 = rsq2, K = K, P = P, power = power,
                  alpha = alpha, test = test, reparameterize = reparameterize) - target)^2
      }
      J <- Jn_optimize(start = min, loss = loss, lower = K + 3, upper = 1e6,
                       solve = "J")
    }
  } else { # solve n
    loss <- function(J) {
      criteria(J = J, n = n, d_est = d_est, d_sd = d_sd, rho_est = rho_est,
               rho_sd = rho_sd, rsq2 = rsq2, K = K, P = P, power = power,
               alpha = alpha, test = test, reparameterize = reparameterize) - target
    }
    min <- 1
    n <- try(stats::uniroot(loss, c(min, 1e8))$root, silent = TRUE)
    # if root-finding method fails, try optimization methods
    if (class(n) == "try-error") {
      loss <- function(n) {
        (criteria(J = J, n = n, d_est = d_est, d_sd = d_sd,
                  rho_est = rho_est, rho_sd = rho_sd, rsq2 = rsq2,
                  K = K, P = P, power = power,
                  alpha = alpha, test = test) - target)^2
      }
      n <- Jn_optimize(start = min, loss = loss, lower = 1, upper = Inf,
                       solve = "n")
    }
  }

  if (plot) {
    ggplot2::theme_set(ggplot2::theme_bw())

    Jn_plots <- plot_Jn(J = J, n = n, d_est = d_est, d_sd = d_sd,
                        rho_est = rho_est, rho_sd = rho_sd,
                        rsq2 = rsq2, K = K, P = P, power = power,
                        alpha = alpha, ep = ep, al = al)
    if (d_sd == 0 & rho_sd == 0) {
      prior_plots <- NULL
    } else {
      prior_plots <- plot_prior(d_est = d_est, d_sd = d_sd,
                                rho_est = rho_est, rho_sd = rho_sd)
    }

    if (J >= 9e5) warning(paste0("Plots may be unreliable."))

    return(list(Jn_plots = Jn_plots, prior_plots = prior_plots,
                Jn = ceiling(cbind(J = J, n = n))))

  } else {
    return(ceiling(cbind(J = J, n = n)))
  }
}
