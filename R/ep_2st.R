#' Expected Power for Independent Sample t-tests
#'
#' \code{ep_2st()} computes the expected power given the effect size estimate
#' and uncertainty in the estimate for an independent sample t-test.
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
#' @return The expected power given the sample size and the prior for the effect size.
#' @export
#' @examples
#' ep_2st(n1 = 100, n2 = 100, d_est = .4, d_sd = .2)
ep_2st <- function(d_est, d_sd, n1, n2, alpha = .05, power = .8,
                   test = "two.sided") {

  # for plotting, assuming n1 = n2
  if (is.null(n2)) n2 <- n1

  if (d_sd == 0) {
    pow_2st(n1 = n1, n2 = n2, d_est = delta, alpha = alpha, test = test)
  } else {
    cubature::hcubature(
      function(delta) {
        pow_2st(n1 = n1, n2 = n2, d_est = delta, alpha = alpha, test = test) *
          stats::dnorm(delta, mean = d_est, sd = d_sd)
      },
      lowerLimit = -Inf, upperLimit = Inf, vectorInterface = TRUE
    )$integral
  }
}
