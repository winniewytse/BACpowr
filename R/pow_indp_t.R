#' Determine Statistical Power for Independent Sample T-Tests
#'
#' \code{pow_indp_t()} computes the statistical power for an independent sample
#'  t-test.
#'
#' @param delta Effect size estimate.
#' @param n1 Sample size of group 1.
#' @param n2 Sample size of group 2.
#' @param alpha Type I error rate. Defaults to \code{.05}.
#' @param test One-tailed or two-tailed test. Defaults to "two.sided".
#' @return Statistical power given the sample size for an independent sample t-test.
#' @export
#' @examples
#' pow_indp_t(n1 = 100, n2 = 100, delta = .4)
pow_indp_t <- function(delta, n1, n2, alpha = .05, test = "two.sided") {
  df <- n1 + n2 - 2
  ncp <- delta * sqrt(n1 * n2 / (n1 + n2))
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
