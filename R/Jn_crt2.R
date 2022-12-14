#' Determine Number of Clusters or Cluster Size for Two-Level CRTs
#'
#' \code{Jn_crt2()} uses the HCB approach to solve for the minimum required
#' number of clusters (J) or cluster size (n) that would achieve a desired
#' expected power or assurance level for a two-level CRT (cluster randomized
#' trial).
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
#'   Defaults to "two-sided".
#' @param plot Whether plots of J and n against expected power or
#'   assurance level should be returned. Defaults to FALSE.
#' @param max_try Maximum number of tries allowed. Defaults to 1e6.
#' @return A 1 x 2 array containing the J and n values that together achieve the
#'   desired expected power or assurance level. If \code{plot = TRUE}, also
#'   returns plots.
#' @import stats methods
#' @export
#' @examples
#' Jn_crt2(delta = .5, delta_sd = .2, rho = .1, rho_sd = .05, J = 30)
#' @seealso \url{https://winnie-wy-tse.shinyapps.io/hcb_shiny/}

Jn_crt2 <- function(delta, delta_sd, rho, rho_sd, rsq2 = 0, J = NULL,
                    n = NULL, K = 0, P = .5, alpha = .05, power = .8,
                    ep = NULL, al = NULL, test = "two.sided",
                    plot = FALSE, max_try = 1e6) {

  if (is.null(J) & is.null(n)) stop(paste0("Please specify either n or J."))

  # If neither EP nor AL is specified, set EP equal to power to solve for power.
  if (is.null(ep) & is.null(al)) ep <- power

  # If both EP and AL were specified, solve sample size for desired AL.
  if (!is.null(ep) & !is.null(al)) ep <- NULL

  params <- list(delta = delta, delta_sd = delta_sd, rho = rho, rho_sd = rho_sd,
                 rsq2 = rsq2, K = K, P = P, power = power, alpha = alpha,
                 test = test)

  Jn_try(J = J, n = n, ep = ep, al = al, params = params, max_try = max_try,
         design = "crt2")

  # As a starting point, compute J and n using the conventional approach.
  Jn_conv <- do.call(Jn_crt2_c, append(list(J = J, n = n),
                                       within(params, rm(delta_sd, rho_sd))))

  # If uncertainty is set to 0 for effect size and ICC estimates, return J and n
  # values computed using the conventional approach and plots if plot == TRUE.
  if (delta_sd == 0 & rho_sd == 0) {
    if (plot) {
      Jn_plots <- do.call(plot_Jn,
                          append(list(J = Jn_conv[1], n = Jn_conv[2]), params))
      return(list(Jn_plots = Jn_plots, Jn = ceiling(Jn_conv)))
    } else {
      return(ceiling(Jn_conv))
    }
  }

  if (!is.null(ep)) { # solve J/n for the target ep
    criteria <- ep_crt2
    target <- ep
    goal <- "ep"
  }
  else if (!is.null(al)) { # solve J/n for the target al
    criteria <- al_crt2
    target <- al
    goal <- "al"
  }

  # Define a loss function for J or n, first attempt with a root-finding method,
  # If root-finding fails, try with optimization methods
  if (is.null(J)) {
    loss_root <- function(J) {
      do.call(criteria, append(list(J = J, n = n), params)) - target
    }
    loss_opt <- function(J) {
      (do.call(criteria, append(list(J = J, n = n), params)) - target)^2
    }
    min <- if (is.null(al) & !is.null(ep)) (K + 2 + 1) else Jn_conv[1]
    given <- "n"
    size <- n
  } else if (is.null(n)) {
    loss_root <- function(n) {
      do.call(criteria, append(list(J = J, n = n), params)) - target
    }
    loss_opt <- loss <- function(n) {
      (do.call(criteria, append(list(J = J, n = n), params)) - target)^2
    }
    min <- 1
    given <- "J"
    size <- J
  }

  message_par <- list(given = given, goal = goal, size = size, target = target)
  root <- try(stats::uniroot(loss_root, interval = c(min, max_try))$root,
              silent = TRUE)
  if(is(root,"try-error")) {
    opt_sol <- optimize_Jn(start = min, loss = loss_opt, lower = min,
                           upper = max_try, message_par = message_par)
    if (is.null(J)) J <- opt_sol
    else if (is.null(n)) n <- opt_sol
  } else {
    if (is.null(J)) J <- root
    else if (is.null(n)) n <- root
  }
  # if (is.null(J)) { # solve for J
  #   # loss <- function(J) {
  #   #   do.call(criteria, append(list(J = J, n = n), params)) - target
  #   # }
  #   min_j <- if (is.null(al) & !is.null(ep)) (K + 2 + 1) else Jn_conv[1]
  #
  #   J <- try(stats::uniroot(loss, interval = c(min_j, 1e8))$root, silent = TRUE)
  #
  #   # if root-finding method fails, try optimization methods
  #   if (class(J) == "try-error") {
  #     loss <- function(J) {
  #       (do.call(criteria, append(list(J = J, n = n), params)) - target)^2
  #     }
  #     J <- optimize_Jn(start = min_j, loss = loss, lower = K + 3, upper = 1e6,
  #                      solve = "J")
  #   }
  # } else if (is.null(n)) { # solve n
  #   # loss <- function(n) {
  #   #   do.call(criteria, append(list(J = J, n = n), params)) - target
  #   # }
  #   min_n <- 1
  #
  #   n <- try(stats::uniroot(loss, interval = c(min_n, 1e8))$root, silent = TRUE)
  #   if (class(n) == "try-error") {
  #     loss <- function(n) {
  #       (do.call(criteria, append(list(J = J, n = n), params)) - target)^2
  #     }
  #     n <- optimize_Jn(start = min_n, loss = loss, lower = 1, upper = Inf,
  #                      solve = "n")
  #   }
  # }

  if (plot) {
    prior_plots <- plot_prior(delta = delta, delta_sd = delta_sd, rho = rho,
                              rho_sd = rho_sd)
    Jn_plots <-  do.call(plot_Jn, append(list(J = J, n = n), params))

    return(list(Jn_plots = Jn_plots, prior_plots = prior_plots,
                Jn = ceiling(cbind(J = J, n = n))))
  } else {
    return(ceiling(cbind(J = J, n = n)))
  }
}
