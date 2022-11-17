#' Determine Assurance Level (AL) for Independent Sample T-Tests
#'
#' \code{al_indp_t()} computes the assurance level of power given the effect size
#'  estimates and the uncertainty level of the parameter estimates for an
#'  independent sample t-test.
#'
#' @param delta Effect size estimate, defined as
#'   \eqn{d = \frac{\bar x_2 - \bar x_1}
#'   {\sqrt{\frac{(n_1 - 1)s_1^2 + (n_2 - 1)s_2^2}{n_1 + n_2 - 2}}}}
#' @param delta_sd Uncertainty level of the effect size estimate.
#' @param n1 Sample size of group 1.
#' @param n2 Sample size of group 2.
#' @param alpha Type I error rate. Defaults to \code{.05}.
#' @param power Desired level of statistical power. Defaults to \code{.8}.
#' @param test One-sided or two-sided test. Defaults to "two.sided".
#' @return Assurance level given the sample size and the prior for the effect size.
#' @export
#' @examples
#' al_indp_t(n1 = 100, n2 = 100, delta = .4, delta_sd = .2)

al_indp_t <- function(n1, n2, delta, delta_sd, alpha = .05, power = .8,
                   test = "two.sided") {

  # for plotting, assuming n1 = n2
  if (is.null(n2)) n2 <- n1

  if (delta_sd == 0) {
    pow_indp_t(n1 = n1, n2 = n2, delta = delta, alpha = alpha, test = test)
  } else {
    d_star <- pow_inv2(power = power, alpha = alpha, n1 = n1, n2 = n2, test = test)
    if (test == "two.sided") {
      stats::pnorm(d_star, mean = delta, sd = delta_sd, lower.tail = FALSE) +
        stats::pnorm(- d_star, mean = delta, sd = delta_sd, lower.tail = TRUE)
    } else if (test == "one.sided") {
      stats::pnorm(d_star, mean = delta, sd = delta_sd, lower.tail = FALSE)
    }
  }
}

pow_inv2 <- function(power, alpha, n1, n2, test) {
  # if (test == "two.sided") alpha_star <- alpha / 2
  # else if (test == "one.sided") alpha_star <- alpha
  # (qnorm(1 - alpha_star) - qnorm(1 - power)) *
  #     sqrt(((n1 - 1) + (n2 - 1)) / (n1 + n2 - 2) * (1 / n1 + 1 / n2))
  df <- n1 + n2 - 2
  if (test == "two.sided") {
    cv <- stats::qt(1 - alpha / 2, df)
    inv <- function(delta) {
      ncp <- delta * sqrt(n1 * n2 / (n1 + n2))
      stats::pt(cv, df, ncp, lower.tail = FALSE) +
        stats::pt(-cv, df, ncp, lower.tail = TRUE) - power
    }
  } else if (test == "one.sided") {
    cv <- stats::qt(1 - alpha, df)
    inv <- function(delta) {
      ncp <- delta * sqrt(n1 * n2 / (n1 + n2))
      stats::pt(cv, df, ncp, lower.tail = FALSE) - power
    }
  }
  stats::uniroot(inv, c(0, 100))$root
}
