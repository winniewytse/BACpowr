#' Determine Number of Clusters or Cluster Size for Two-Level MSRTs
#'
#' \code{Jn_msrt2()} uses the HCB approach to solve for the minimum required
#' number of clusters (J) or cluster size (n) that would achieve a desired
#' expected power or assurance level for a two-level MSRT (multisite randomized
#' control trial).
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
#'    \eqn{\sigma^2} are the variance components in the unconditional model.
#' @param rho_sd Uncertainty level of the ICC estimate.
#' @param omega Estimate of the treatment effect heterogeneity, defined as
#'   \eqn{\omega = \frac{\tau_1^2}{\tau_0^2}} where \eqn{\tau_0^2} is the
#'   variance of the intercept random component and \eqn{\tau_1^2} is the
#'   variance of the treatment random effect.
#' @param omega_sd Uncertainty level of the treatment effect heterogeneity
#'   estimate.
#' @param rsq1 Estimate of variance explained by the level-1 (individual-level)
#'   covariates.
#' @param rsq2 Estimate of variance explained by the level-2 (cluster-level)
#'  covariates.
#' @param J Number of clusters. If \code{J} is specified, \code{Jn_msrt2()}
#'   determines \code{n} for the specified expected power/assurance level.
#' @param n Cluster size. If \code{n} is specified, \code{Jn_msrt2()}
#'   determines \code{J} for the specified expected power/assurance level.
#' @param K Number of cluster-level covariates. \code{0} by default.
#' @param P Proportion of the clusters that is treatment group. \code{.5} by
#'   default.
#' @param power Desired level of statistical power. \code{.8} by default.
#' @param alpha Type I error rate. \code{.05} by default.
#' @param ep Desired expected power. If neither \code{ep} nor \code{al} is
#'   specified, \code{ep} = \code{power}. For example, an 80% EP indicates that
#'   the mean or average power value is 80% over the specified uncertainty.
#' @param al Desired assurance level. For example, an 80% AL indicates 80%
#'   of the power values are above the desired statistical power over the
#'   specified uncertainty.
#' @param test Whether a one-sided or two-sided test should be performed.
#'   Defaults to "two-sided".
#' @param plot Whether plots of J and n against expected power or
#'   assurance level should be returned. Defaults to FALSE.
#' @return A 1 x 2 array containing the J and n values that together achieve the
#'   desired expected power or assurance level. If \code{plot = TRUE}, also
#'   returns plots.
#' @import stats
#' @export
#' @examples
#' Jn_msrt2(delta = .5, delta_sd = .1, rho = .1, rho_sd = .1, omega = .3,
#'          omega_sd = .1, J = 30)
#' Jn_msrt2(delta = .5, delta_sd = 0, rho = .1, rho_sd = 0, omega = .3,
#'          omega_sd = 0, n = 5)
#' @seealso \url{https://winnie-wy-tse.shinyapps.io/hcb_shiny/}

Jn_msrt2 <- function(delta, delta_sd, rho, rho_sd, omega, omega_sd, rsq1 = 0,
                     rsq2 = 0, J = NULL, n = NULL, K = 0, P = .5, alpha = .05,
                     power = .8, ep = NULL, al = NULL, test = "two.sided",
                     plot = FALSE, max_try = 1e6) {

  if (is.null(J) & is.null(n)) stop(paste0("Please specify either n or J."))

  # If neither EP nor AL is specified, set EP equal to power to solve for power.
  if (is.null(ep) & is.null(al)) ep <- power

  params <- list(delta = delta, delta_sd = delta_sd, rho = rho, rho_sd = rho_sd,
                 omega = omega, omega_sd = omega_sd, rsq1 = rsq1, rsq2 = rsq2,
                 K = K, P = P, power = power, alpha = alpha, test = test)

  Jn_try(J = J, n = n, ep = ep, al = al, params = params, max_try = max_try,
         design = "msrt2")

  # use Jn with the conventional approach as starting points for efficiency
  sol_c <- Jn_msrt2_c(delta = delta, rho = rho, omega = omega, rsq1 = rsq1,
                      rsq2 = rsq2, J = J, n = n, K = K, P = P, alpha = alpha,
                      power = power, test = test)
  if (delta_sd == 0 & rho_sd == 0 & omega_sd == 0) {
    if (plot) {
      Jn_plots <- plot_Jn(J = sol_c[1], n = sol_c[2],
                          delta = delta, delta_sd = delta_sd, rho = rho,
                          rho_sd = rho_sd, omega = omega, omega_sd = omega_sd,
                          rsq1 = rsq1, rsq2 = rsq2, K = K, P = P, power = power,
                          alpha = alpha, ep = ep, al = al)
      return(list(Jn_plots = Jn_plots, Jn = ceiling(sol_c)))
    } else {
      return(ceiling(sol_c))
    }
  }

  if (is.null(al)) { # solve with the expected power
    criteria <- ep_msrt2
    target <- ep
    goal <- "ep"
  } else { # solve with the assurance level
    criteria <- al_msrt2
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
    min <- K + 1 + 1
    start <- sol_c[1]
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
    start <- sol_c[2]
    given <- "J"
    size <- J
  }
  message_par <- list(given = given, goal = goal, size = size, target = target)
  root <- try(stats::uniroot(loss_root, interval = c(min, max_try))$root,
              silent = TRUE)
  if (class(root) == "try-error") {
    opt_sol <- optimize_Jn(start = start, loss = loss_opt, lower = min,
                           upper = max_try, message_par = message_par)
    if (is.null(J)) J <- opt_sol
    else if (is.null(n)) n <- opt_sol
  } else {
    if (is.null(J)) J <- root
    else if (is.null(n)) n <- root
  }

  if (J >= 9e5) warning(paste0("Results may be unreliable due to convergence
                                 issues. Try increasing the cluster size (n),
                                 using smaller uncertainty (delta_sd or rho_sd),
                                 or decreasing EP or AL."))
  if (plot) {
    if (sum(c(delta_sd, rho_sd, omega_sd) != 0)) smooth <- 21 else smooth <- 51

    Jn_plots <- do.call(plot_Jn, append(list(J = J, n = n, smooth = smooth),
                                        params))

    prior_plots <- plot_prior(delta = delta, delta_sd = delta_sd, rho = rho,
                              rho_sd = rho_sd, omega = omega,
                              omega_sd = omega_sd)

    return(list(Jn_plots = Jn_plots, prior_plots = prior_plots,
                Jn = ceiling(cbind(J = J, n = n))))
  } else {
    return(ceiling(cbind(J = J, n = n)))
  }
}
