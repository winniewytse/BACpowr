#' Assurance Level for Independent Sample t-tests.
#'
#' \code{al_crt()} computes the assurance level of power given the estimates and
#' the uncertainty level of the parameter estimates for an independent sample t-test.
#'
#' @param d_est Effect size estimate, defined as
#'   \eqn{d = \frac{\bar x_2 - \bar x_1}
#'   {\sqrt{\frac{(n_1 - 1)s_1^2 + (n_2 - 1)s_2^2}{n_1 + n_2 - 2}}}}
#' @param d_sd Uncertainty level of the effect size estimate.
#' @param n1 Sample size of group 1.
#' @param n2 Sample size of group 2.
#' @param alpha Type I error rate. Default to be \code{.05}.
#' @param power Desired statistical power to achieve. Default to be \code{.8}.
#' @param test One-sided or two-sided test. Options are either "one.sided" or "two.sided".
#' @return The assurance level given the sample size and the prior for the effect size.
#' @export
#' @examples
#' al_2st(n1 = 100, n2 = 100, d_est = .4, d_sd = .2)

# define assurance level function
al_2st <- function(n1, n2, d_est, d_sd, alpha = .05, power = .8,
                   prior_d = c("norm", "trunc_norm"), trunc_d = c(-Inf, Inf),
                   test = "two.sided") {

  # for plotting, assuming n1 = n2
  if (is.null(n2)) n2 <- n1

  if (d_sd == 0) {
    pow_2st(n1 = n1, n2 = n2, d_est = d_est, alpha = alpha, test = test)
  } else {
    if (prior_d == "norm") {
      if (test == "two.sided") {
        d_star <- pow_inv2(power = power, alpha = alpha,
                           n1 = n1, n2 = n2, test = test)
        stats::pnorm(d_star, mean = d_est, sd = d_sd, lower.tail = FALSE) +
          stats::pnorm(- d_star, mean = d_est, sd = d_sd, lower.tail = TRUE)
      } else if (test == "one.sided") {
        stats::pnorm(d_star, mean = d_est, sd = d_sd, lower.tail = FALSE)
      }
    } else if (prior_d == "trunc_norm") {
      d_star <- pow_inv2_trunc(power = power, alpha = alpha, trunc_d = trunc_d,
                               n1 = n1, n2 = n2, test = test)
      if (test == "two.sided") {
        1 - (ptruncnorm(d_star, a = trunc_d[1], b = trunc_d[2],
                        mean = d_est, sd = d_sd, lower.tail = FALSE) +
               ptruncnorm(- d_star, a = trunc_d[1], b = trunc_d[2],
                          mean = d_est, sd = d_sd, lower.tail = TRUE))
      } else if (test == "one.sided") {
        1 - (ptruncnorm(d_star, a = trunc_d[1], b = trunc_d[2],
                        mean = d_est, sd = d_sd, lower.tail = FALSE))
      }
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
    inv <- function(d_est) {
      ncp <- d_est * sqrt(n1 * n2 / (n1 + n2))
      stats::pt(cv, df, ncp, lower.tail = FALSE) +
        stats::pt(-cv, df, ncp, lower.tail = TRUE) - power
    }
  } else if (test == "one.sided") {
    cv <- stats::qt(1 - alpha, df)
    inv <- function(d_est) {
      ncp <- d_est * sqrt(n1 * n2 / (n1 + n2))
      stats::pt(cv, df, ncp, lower.tail = FALSE) - power
    }
  }
  stats::uniroot(inv, c(0, 100))$root
}

pow_inv2_trunc <- function(power, alpha, n1, n2, test, trunc_d) {
  df <- n1 + n2 - 2
  # trunc_t <- d2t_2st(trunc_d, n1, n2, s1 = 1, s2 = 1)
  trunc_ncp <- trunc_d * sqrt(n1 * n2 / (n1 + n2))
  if (test == "two.sided") {
    cv <- stats::qt(1 - alpha / 2, df)
    # cv <- qtrunct(1 - alpha / 2, a = trunc_t[1], b = trunc_t[2], df = df)
    inv <- function(ncp) {
      stats::pt(cv, df, ncp, lower.tail = FALSE) +
        stats::pt(-cv, df, ncp, lower.tail = TRUE) - power
    }
  } else if (test == "one.sided") {
    cv <- stats::qt(1 - alpha / 2, df)
    inv <- function(ncp) {
      stats::pt(cv, df, ncp, lower.tail = FALSE) - power
    }
  }
  ncp_sol <- stats::uniroot(inv, c(trunc_ncp[1], 100))$root
  ncp_sol / sqrt(n1 * n2 / (n1 + n2))
}

d2t_2st <- function(d, n1, n2, s1, s2) {
  d * sqrt(((n1 - 1) * s1^2 + (n2 - 1) * s2^2) / (n1 + n2 - 2)) *
    sqrt(n1 * n2 / (n2 * s1^2 + n1 * s2^2))
}
