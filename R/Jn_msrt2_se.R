#' Determine Number of Clusters or Cluster Size for Two-Level Multisite Randomized Trials
#'
#' @param rho Intraclass correlation value, defined as
#'   \eqn{\rho = \frac{\tau^2}{\tau^2 + \sigma^2}}, where \eqn{\tau^2} and \eqn{\sigma^2}
#'   are the variance components in the unconditional model.
#' @param rho_sd Uncertainty level of the intraclass correlation value.
#' @param omega Treatment effect hetereogeneity value, defined as
#'   \eqn{\omega = \frac{\tau_1^2}{\tau_0^2}} where \eqn{\tau_0^2} is the variance of the
#'   intercept random component and \eqn{\tau_1^2} is the variance of the treatment
#'   random effect.
#' @param omega_sd Uncertainty level of the treatment effect hetereogeneity value.
#' @param rsq1 Estimate of variance explained by the level-1 (e.g., individual-level) covariates.
#' @param rsq2 Estimate of variance explained by the cluster-level covariates.
#' @param J Number of clusters. Determine \code{n} if \code{J} is specified.
#' @param n Cluster size. Determine \code{J} if \code{n} is specified.
#' @param K Number of cluster-level covariates.
#' @param P Proportion of the clusters that is treatment group.
#' @param ese Desired expected confidence interval half width to achieve.
#'   If neither \code{ese} nor \code{ase} is specified, \code{ese} = \code{se}.
#'   An 0.05 expected half width indicates that the average confidence interval
#'   half width is 0.05 over the specified uncertainty.
#' @param ase Assurance level of confidence interval half width to achieve.
#'   An 80% assurance level indicates there is an 80% chance that the confidence interval
#'   half width is below or at the desired half width over the specified uncertainty.
#' @param plot If TRUE, plots of J and n against the expected power or assurance level
#'   will be printed.
#' @return The required J or n and optional plots that show the curves of
#'   expected power/assurance level.
#' @import stats
#' @export
#' @examples
#' Jn_msrt2_se(rho = .1, rho_sd = .1, omega = .3, omega_sd = .1, J = 30, se = .05)
#' Jn_msrt2_se(rho = .1, rho_sd = 0, omega = .3, omega_sd = 0, n = 5, se = .05, ase = .6)
#' @seealso \url{https://winnie-wy-tse.shinyapps.io/hcb_shiny/}

Jn_msrt2_se <- function(rho, rho_sd, omega, omega_sd,
                        rsq1 = 0, rsq2 = 0, J = NULL, n = NULL, K = 0, P = .5,
                        se = .05, ese = NULL, ase = NULL, plot = FALSE) {

  ggplot2::theme_set(ggplot2::theme_bw())

  if (is.null(ese) & is.null(ase)) ese <- se

  # use Jn with the conventional approach as starting points for efficiency
  Jn_msrt <- Jn_msrt2_se_c(rho = rho, omega = omega,
                           rsq1 = rsq1, rsq2 = rsq2,
                           J = J, n = n, K = K, P = P, se = se)
  if (rho_sd == 0 & omega_sd == 0) {
    if (plot) {
      warning("Plots not supported yet. ")
      # Jn_plots <- plot_Jn(J = Jn_msrt[1], n = Jn_msrt[2],
      #                     d_est = d_est, d_sd = d_sd,
      #                     rho_est = rho_est, rho_sd = rho_sd,
      #                     omega_est = omega_est, omega_sd = omega_sd,
      #                     rsq1 = rsq1, rsq2 = rsq2, K = K, P = P,
      #                     power = power, alpha = alpha, ep = ep, al = al)
      # return(list(Jn_plots = Jn_plots, Jn = round(Jn_msrt)))
      return(round(Jn_msrt))
    } else {
      return(round(Jn_msrt))
    }
  }

  if (is.null(ase)) { # solve with the expected power
    criteria <- ese_msrt2
    target <- ese
    if (is.null(J)) {
      min <- K + 2 + 1
    } else {
      min <- Jn_msrt[2]
    }
  } else { # solve with the assurance level
    criteria <- ase_msrt2
    target <- ase
    if (is.null(J)) {
      # set a higher min J to avoid being stuck at the local minimum
      min <- Jn_msrt[1]
    } else {
      min <- Jn_msrt[2]
    }
  }

  if (is.null(J)) { # solve J
    loss <- function(J) {
      criteria(J = J, n = n, rho = rho, rho_sd = rho_sd,
               omega = omega, omega_sd = omega_sd,
               rsq1 = rsq1, rsq2 = rsq2, se = se, K = K, P = P) - target
    }
    J <- try(stats::uniroot(loss, c(min, 1e8))$root, silent = TRUE)
    # if root-finding method fails, try optimization methods
    if (class(J) == "try-error") {
      loss <- function(J) {
        (criteria(J = J, n = n, rho = rho, rho_sd = rho_sd,
                  omega = omega, omega_sd = omega_sd,
                  rsq1 = rsq1, rsq2 = rsq2, se = se, K = K, P = P) - target)^2
      }
      J <- Jn_optimize(start = min, loss = loss, lower = K + 3, upper = 1e6,
                       solve = "J")
    }
  } else { # solve n
    loss <- function(n) {
      criteria(J = J, n = n, rho = rho, rho_sd = rho_sd,
               omega = omega, omega_sd = omega_sd,
               rsq1 = rsq1, rsq2 = rsq2, se = se, K = K, P = P) - target
    }
    n <- try(stats::uniroot(loss, c(min, 1e8))$root, silent = TRUE)
    # if root-finding method fails, try optimization methods
    if (class(n) == "try-error") {
      loss <- function(n) {
        (
          criteria(J = J, n = n, rho = rho, rho_sd = rho_sd,
                  omega = omega, omega_sd = omega_sd,
                  rsq1 = rsq1, rsq2 = rsq2, se = se, K = K, P = P) - target)^2
      }
      n <- Jn_optimize(start = min, loss = loss, lower = 1, upper = Inf,
                       solve = "n")
    }
  }

  if (plot) {
    warning("Plots not supported yet. ")
    # if (sum(c(d_sd, rho_sd, omega_sd) != 0)) smooth <- 21
    # else smooth <- 51
    # Jn_plots <- plot_Jn(J = J, n = n, d_est = d_est, d_sd = d_sd,
    #                     rho_est = rho_est, rho_sd = rho_sd,
    #                     omega_est = omega_est, omega_sd = omega_sd,
    #                     rsq1 = rsq1, rsq2 = rsq2,
    #                     K = K, P = P, power = power,alpha = alpha,
    #                     ep = ep, al = al, smooth = smooth)
    # if (d_sd == 0 & rho_sd == 0 & omega_sd == 0) {
    #   prior_plots <- NULL
    # } else {
    #   prior_plots <- plot_prior(d_est = d_est, d_sd = d_sd,
    #                             rho_est = rho_est, rho_sd = rho_sd,
    #                             omega_est = omega_est, omega_sd = omega_sd)
    # }
    #
    # if (J >= 9e5) warning(paste0("Plots may be unreliable."))
    #
    # return(list(Jn_plots = Jn_plots, prior_plots = prior_plots,
    #             Jn = round(cbind(J = J, n = n))))
    return(round(cbind(J = J, n = n)))

  } else {
    return(round(cbind(J = J, n = n)))
  }
}

# Solve Jn using the conventional approach
Jn_msrt2_se_c <- function(rho, omega, rsq1 = 0, rsq2 = 0,
                          J = NULL, n = NULL, K = 0, P = .5,
                          alpha = .05, se = .05) {

  if (is.null(J)) { # solve J
    loss <- function(J) {
      se_msrt2(J = J, n = n, rho = rho, omega = omega,
               rsq1 = rsq1, rsq2 = rsq2, K = K, P = P) - se
    }
    min <- K + 2 + 1
    J <- try(stats::uniroot(loss, c(min, 1e8))$root, silent = TRUE)
    if (class(J) == "try-error") {
      loss <- function(J) {
        (se_msrt2(J = J, n = n, rho = rho, omega = omega,
                  rsq1 = rsq1, rsq2 = rsq2, K = K, P = P) - se)^2
      }
      J <- try(Jn_optimize(start = min, loss = loss, lower = K + 3, upper = 1e6),
               silent = TRUE)
      if (class(J) == "try-error") {
        stop(paste0("The specified cluster size (n) is insufficient to achieve ",
                    "the desired half width (standard error) of the effect size. ",
                    "Please increase n. "))
      }
    }
  } else { # solve n
    loss <- function(n) {
      se_msrt2(J = J, n = n, rho = rho, omega = omega,
               rsq1 = rsq1, rsq2 = rsq2, K = K, P = P) - se
    }
    min <- 1
    n <- try(stats::uniroot(loss, c(min, 1e8))$root, silent = TRUE)
    if (class(n) == "try-error") {
      stop(paste0("The specified number of clusters (J) is insufficient to achieve ",
                  "the desired half width (standard error) of the effect size. ",
                  "Please increase J. "))
    }
  }
  return(cbind(J = J, n = n))
}
