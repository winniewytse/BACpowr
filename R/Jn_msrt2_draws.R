#' Determine Number of Clusters or Cluster Size for Two-Level Multisite Randomzied Trials
#'
#' @export

Jn_msrt2_draws <- function(draws_d, draws_rho, draws_omega,
                           rsq1 = 0, rsq2 = 0,
                           J = NULL, n = NULL, K = 0, P = .5,
                           alpha = .05, power = .8, precision = .15,
                           ep = NULL, al = NULL, apr = NULL,
                           test = "two.sided", plot = FALSE) {

  ggplot2::theme_set(ggplot2::theme_bw())

  # use Jn with the conventional approach as starting points for efficiency
  Jn_msrt <- Jn_msrt2_c(delta = mean(draws_d), rho = get_mode(draws_rho),
                        omega = get_mode(draws_omega),
                        rsq1 = rsq1, rsq2 = rsq2, J = J, n = n, K = K, P = P,
                        alpha = alpha, power = power, test = test)

  if (!is.null(al)) { # solve with the assurance level
    goal <- "al"
    target <- al
    if (is.null(J)) {
      # set a higher min J to avoid being stuck at the local minimum
      min <- Jn_msrt[1]
    }
  } else if (!is.null(ep)) { # solve with the expected power
    goal <- "ep"
    target <- ep
    if (is.null(J)) {
      min <- K + 2 + 1
    }
  } else if (!is.null(apr)) {
    goal = "apr"
    target <- apr
    if (is.null(J)) {
      min <- K + 2 + 1
    }
  }

  if (is.null(J)) { # solve J
    loss <- function(J) {
      draws_msrt2(draws_d = draws_d, draws_rho = draws_rho,
                  draws_omega = draws_omega, rsq1 = rsq1, rsq2 = rsq2,
                  J = J, n = n, K = K, P = P,
                  power = power, precision = precision,
                  goal = goal, alpha = alpha, test = test) - target
    }
    J <- try(stats::uniroot(loss, c(min, 1e8))$root, silent = TRUE)
    # if root-finding method fails, try optimization methods
    if (class(J) == "try-error") {
      loss <- function(J) {
        (draws_msrt2(draws_d = draws_d, draws_rho = draws_rho,
                     draws_omega = draws_omega, rsq1 = rsq1, rsq2 = rsq2,
                     J = J, n = n, K = K, P = P,
                     power = power, precision = precision,
                     goal = goal, alpha = alpha, test = test) - target)^2
      }
      J <- Jn_optimize(start = min, loss = loss, lower = K + 3, upper = 1e6,
                       solve = "J")
    }
  } else { # solve n
    loss <- function(n) {
      draws_msrt2(draws_d = draws_d, draws_rho = draws_rho,
                  draws_omega = draws_omega, rsq1 = rsq1, rsq2 = rsq2,
                  J = J, n = n, K = K, P = P,
                  power = power, precision = precision,
                  goal = goal, alpha = alpha, test = test) - target
    }
    min <- 1
    n <- try(stats::uniroot(loss, c(min, 1e8))$root, silent = TRUE)
    # if root-finding method fails, try optimization methods
    if (class(n) == "try-error") {
      loss <- function(n) {
        (draws_msrt2(draws_d = draws_d, draws_rho = draws_rho,
                     draws_omega = draws_omega, rsq1 = rsq1, rsq2 = rsq2,
                     J = J, n = n, K = K, P = P,
                     power = power, precision = precision,
                     goal = goal, alpha = alpha, test = test) - target)^2
      }
      n <- Jn_optimize(start = min, loss = loss, lower = 1, upper = Inf,
                       solve = "n")
    }
  }

  # if (plot) {
  #   if (sum(c(delta_sd, rho_sd, omega_sd) != 0)) smooth <- 21
  #   else smooth <- 51
  #   Jn_plots <- plot_Jn(J = J, n = n, delta = delta, delta_sd = delta_sd,
  #                       rho = rho, rho_sd = rho_sd,
  #                       omega = omega, omega_sd = omega_sd,
  #                       rsq1 = rsq1, rsq2 = rsq2,
  #                       K = K, P = P, power = power,alpha = alpha,
  #                       ep = ep, al = al, smooth = smooth)
  #   if (delta_sd == 0 & rho_sd == 0 & omega_sd == 0) {
  #     prior_plots <- NULL
  #   } else {
  #     prior_plots <- plot_prior(delta = delta, delta_sd = delta_sd,
  #                               rho = rho, rho_sd = rho_sd,
  #                               omega = omega, omega_sd = omega_sd)
  #   }
  #
  #   if (J >= 9e5) warning(paste0("Plots may be unreliable."))
  #
  #   return(list(Jn_plots = Jn_plots, prior_plots = prior_plots,
  #               Jn = ceiling(cbind(J = J, n = n))))
  #
  # } else {
  return(ceiling(cbind(J = J, n = n)))
  # }
}


get_mode <- function(draws) {
  dens <- density(draws)
  dens$x[which.max(dens$y)]
}

