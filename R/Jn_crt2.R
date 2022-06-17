#' Determine Number of Clusters or Cluster Size for Two-Level CRTs
#'
#' \code{Jn_crt2()} solves for the required number of clusters (J) or cluster size (n)
#' for a two-level CRT using the HCB approach. When the uncertainty level of a parameter
#' is specified, this function determines the minimum required sample size that achieves
#' the desired expected power or assurance level.
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
#' @param rsq2 Estimate of variance explained by the cluster-level covariates.
#' @param J Number of clusters. Determine \code{n} if \code{J} is specified.
#' @param n Cluster size. Determine \code{J} if \code{n} is specified.
#' @param K Number of cluster-level covariates.
#' @param P Proportion of the clusters that is treatment group.
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
#' @return The required J or n and optional plots that show the curves of
#'   expected power/assurance level.
#' @import stats
#' @export
#' @examples
#' Jn_crt2(d_est = .5, d_sd = .2, rho_est = .1, rho_sd = .05, J = 30)
#' @seealso \url{https://winnie-wy-tse.shinyapps.io/hcb_shiny/}

Jn_crt2 <- function(d_est, d_sd, rho_est, rho_sd, rsq2 = 0,
                    J = NULL, n = NULL, K = 0, P = .5,
                    alpha = .05, power = .8, ep = NULL, al = NULL,
                    test = "two.sided", plot = FALSE){

  if (is.null(ep) & is.null(al)) ep <- power

  # use Jn with the conventional approach as starting points for efficiency
  Jn_crt <- Jn_crt2_c(d_est = d_est, rho_est = rho_est, rsq2 = rsq2,
                      J = J, n = n, K = K, P = P,
                      alpha = alpha, power = power, test = test)
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

  if (is.null(al)) { # solve with the expected power
    criteria <- ep_crt2
    target <- ep
    if (is.null(J)) {
      min <- K + 2 + 1
    }
  } else { # solve with the assurance level
    criteria <- al_crt2
    target <- al
    if (is.null(J)) {
      # set a higher min J to avoid being stuck at the local minimum
      min <- Jn_crt[1]
    }
  }

  if (is.null(J)) { # solve J
    loss <- function(J) {
      criteria(J = J, n = n, d_est = d_est, d_sd = d_sd,
               rho_est = rho_est, rho_sd = rho_sd, rsq2 = rsq2,
               K = K, P = P, power = power, alpha = alpha, test = test) - target
    }
    J <- try(stats::uniroot(loss, c(min, 1e8))$root, silent = TRUE)
    # if root-finding method fails, try optimization methods
    if (class(J) == "try-error") {
      loss <- function(J) {
        (criteria(J = J, n = n, d_est = d_est, d_sd = d_sd,
                  rho_est = rho_est, rho_sd = rho_sd, rsq2 = rsq2,
                  K = K, P = P, power = power,
                  alpha = alpha, test = test) - target)^2
      }
      J <- Jn_optimize(start = min, loss = loss, lower = K + 3, upper = 1e6,
                       solve = "J")
    }
  } else { # solve n
    loss <- function(n) {
      criteria(J = J, n = n, d_est = d_est, d_sd = d_sd,
               rho_est = rho_est, rho_sd = rho_sd, rsq2 = rsq2,
               K = K, P = P, power = power, alpha = alpha, test = test) - target
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
    prior_plots <- plot_prior(d_est = d_est, d_sd = d_sd,
                              rho_est = rho_est, rho_sd = rho_sd)

    if (J >= 9e5) warning(paste0("Plots may be unreliable."))

    return(list(Jn_plots = Jn_plots, prior_plots = prior_plots,
                Jn = ceiling(cbind(J = J, n = n))))

  } else {
    return(ceiling(cbind(J = J, n = n)))
  }
}

# Solve J/n by optimization method
Jn_optimize <- function(start, loss, lower, upper, solve) {
  # try using PORT routines
  port <- stats::nlminb(start = start, objective = loss, lower = lower)
  # if PORT routines do not converge, try Brent
  if (port$objective > 1e-5) {
    brent <- try(stats::optim(par = start, fn = loss, lower = lower,
                              upper = upper, method = "Brent"),
                 silent = TRUE)
    # if Brent does not work, try LBFGSB
    if (class(brent) == "try-error") {
      condition <- "error"
    } else {
      if (brent$value > 1e-5) {
        condition <- "non-convergence"
      } else {
        condition <- "convergence"
      }
    }
    if (condition %in% c("error", "non-convergence")) {
      lbfgsb <- stats::optim(par = start, fn = loss, lower = lower,
                             upper = Inf, method = "L-BFGS-B")
      if (lbfgsb$value > 1e-5) {
        sol <- lbfgsb$par
        if (solve == "n") {
          warning("The algorithm fails to converge for the specified priors. ",
                  "Please consider increasing J or reducing the expected power or",
                  "assurance level. ")
        } else {
          warning(paste0("The algorithm fails to converge for the specified priors. ",
                         "There may not exist a solution for the desired expected ",
                         "power or assurance level. ",
                         "Please consider some lower power/assurance level. "))
        }
      } else {
        sol <- lbfgsb$par
      }
    } else {
      sol <- brent$par
    }
  } else {
    sol <- port$par
  }
  return(sol)
}

# Solve Jn using the conventional approach
Jn_crt2_c <- function(d_est, rho_est, rsq2 = 0,
                      J = NULL, n = NULL, K = 0, P = .5,
                      alpha = .05, power = .80, test = "two.sided") {

  if (is.null(J)) { # solve J
    loss <- function(J) {
      pow_crt2(J = J, n = n, d_est = d_est, rho_est = rho_est,
               rsq2 = rsq2, test = test, P = P) - power
    }
    min <- K + 2 + 1
    J <- try(stats::uniroot(loss, c(min, 1e8))$root, silent = TRUE)
    if (class(J) == "try-error") {
      loss <- function(J) {
        (pow_crt2(J = J, n = n, d_est = d_est, rho_est = rho_est,
                  rsq2 = rsq2, test = test, P = P) - power)^2
      }
      J <- Jn_optimize(start = min, loss = loss, lower = K + 3, upper = 1e6)
    }
  } else { # solve n
    loss <- function(n) {
      pow_crt2(J = J, n = n, d_est = d_est, rho_est = rho_est,
               rsq2 = rsq2, test = test, P = P) - power
    }
    min <- 1
    n <- stats::uniroot(loss, c(min, 1e8))$root
  }
  return(cbind(J = J, n = n))
}
