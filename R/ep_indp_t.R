#' Determine Expected Power (EP) for Independent Sample T-Tests
#'
#' \code{ep_indp_t()} computes the expected power given the effect size estimate
#' and uncertainty in the estimate for an independent sample t-test.
#'
#' @param delta Effect size estimate, defined as
#'   \eqn{d = \frac{\bar x_2 - \bar x_1}
#'   {\sqrt{\frac{(n_1 - 1)s_1^2 + (n_2 - 1)s_2^2}{n_1 + n_2 - 2}}}}
#' @param delta_sd Uncertainty level of the effect size estimate.
#' @param n1 Sample size of group 1.
#' @param n2 Sample size of group 2.
#' @param alpha Type I error rate. Defaults to \code{.05}.
#' @param power Desired statistical power to achieve. Defaults to \code{.8}.
#' @param test One-sided or two-sided test. Defaults to "two.sided".
#' @return Expected power given the sample size and the prior for the effect size.
#' @export
#' @examples
#' ep_indp_t(n1 = 100, n2 = 100, delta = .4, delta_sd = .2)
ep_indp_t <- function(delta, delta_sd, n1, n2, alpha = .05, power = .8,
                   test = "two.sided") {

  # for plotting, assuming n1 = n2
  if (is.null(n2)) n2 <- n1

  if (delta_sd == 0) {
    pow_indp_t(n1 = n1, n2 = n2, delta = delta, alpha = alpha, test = test)
  } else {
    cubature::hcubature(
      function(delta) {
        pow_indp_t(n1 = n1, n2 = n2, delta = delta, alpha = alpha, test = test) *
          stats::dnorm(delta, mean = delta, sd = delta_sd)
      },
      lowerLimit = -Inf, upperLimit = Inf, vectorInterface = TRUE
    )$integral
  }
}
