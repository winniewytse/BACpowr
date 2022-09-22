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
                   prior_d = c("norm", "trunc_norm", "zero-inflated"),
                   trunc_d = c(-Inf, Inf),
                   test = "two.sided", ndraws = 1e6) {

  # for plotting, assuming n1 = n2
  if (is.null(n2)) n2 <- n1

  if (d_sd == 0) {
    pow_2st(n1 = n1, n2 = n2, d_est = d_est, alpha = alpha, test = test)
  } else {
    if (prior_d == "zero-inflated") {
      df <- n1 + n2 - 2
      cv <- stats::qt(1 - alpha / 2, df)
      d_draws <- rnorm(ndraws, mean = d_est, sd = d_sd)
      d_draws[d_draws < 0] <- 0
      ncp_draws <- d_draws * sqrt(n1 * n2 / (n1 + n2))
      pow_draws <- stats::pt(cv, df = df, ncp = ncp_draws, lower.tail = FALSE) +
        stats::pt(-cv, df = df, ncp = ncp_draws, lower.tail = TRUE)
      mean(pow_draws)
    } else if (prior_d == "norm") {
      cubature::hcubature(
        function(delta) {
          pow_2st(n1 = n1, n2 = n2, d_est = delta, alpha = alpha, test = test) *
            stats::dnorm(delta, mean = d_est, sd = d_sd)
        },
        lowerLimit = -Inf, upperLimit = Inf, vectorInterface = TRUE
      )$integral
    } else if (prior_d == "trunc_norm") {
      cubature::cuhre(
        function(delta) {
          pow_2st(n1 = n1, n2 = n2, d_est = delta, alpha = alpha, test = test) *
            dtruncnorm(delta, a = trunc_d[1], b = trunc_d[2],
                       mean = d_est, sd = d_sd)
        },
        lowerLimit = trunc_d[1], upperLimit = trunc_d[2],
      )$integral
    }
  }
}
