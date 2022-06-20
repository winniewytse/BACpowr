#' Statistical Power for Independent Sample t-tests
#'
#' \code{pow_2st()} computes the statistical power for an independent sample t-test.
#'
#' @param d_est Effect size estimate.
#' @param n1 Sample size of group 1.
#' @param n2 Sample size of group 2.
#' @param alpha Type I error rate. Default to be \code{.05}.
#' @param test One-tailed or two-tailed test.
#' @return Statistical power given the sample size for an independent sample t-test.
#' @export
#' @examples
#' pow_2st(n1 = 100, n2 = 100, d_est = .4)
pow_2st <- function(d_est, n1, n2, alpha = .05, test = "two.sided") {
  df <- n1 + n2 - 2
  ncp <- d_est * sqrt(n1 * n2 / (n1 + n2))
  if (test == "two.sided") {
    cv <- stats::qt(1 - alpha / 2, df)
    pow <- stats::pt(cv, df = df, ncp = ncp, lower.tail = FALSE) +
      stats::pt(-cv, df = df, ncp = ncp, lower.tail = TRUE)
  } else if (test == "one.sided") {
    cv <- stats::qt(1 - alpha, df)
    pow <- stats::pt(cv, df = df, ncp = ncp, lower.tail = FALSE)
  }
  return(pow)
}
