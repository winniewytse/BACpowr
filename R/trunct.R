#' Truncated t distribution
#' @export

dtrunct <- function(x, a = -Inf, b = Inf, df, ncp = 0, ...) {
  if (x <= a | x >= b) {
    0
  } else {
    dt(x, df = df, ncp = ncp, ...) /
      (pt(b, df = df, ncp = ncp, ...) - pt(a, df = df, ncp = ncp, ...))
  }
}
# cdf of truncated normal
ptrunct <- function(q, a = -Inf, b = Inf, df, ncp = 0, ...) {
  if (q <= a) 0
  else if (q >= b) 1
  else
    (pt(q, df = df, ncp = ncp, ...) - pt(a, df = df, ncp = ncp, ...)) /
    (pt(b, df = df, ncp = ncp, ...) - pt(a, df = df, ncp = ncp, ...))
}
# quantile function of truncated normal
qtrunct <- function(p, a = -Inf, b = Inf, df, ncp = 0, ...) {
  qt(pt(a, df = df, ncp = ncp, ...) +
       p * (pt(b, df = df, ncp = ncp, ...) - pt(a, df = df, ncp = ncp, ...)),
     df = df, ncp = ncp, ...)
}
# random draws of truncated normal
rtrunct <- function(n, a = -Inf, b = Inf, df, ncp = 0) {
  ndraws <- 0
  draws_ab <- c()
  while (ndraws < n) {
    draws <- rt(n, df = df, ncp = ncp, ...)
    draws_ab <- c(draws_ab, draws[draws > a & draws < b])
    ndraws <- length(draws_ab)
  }
  return(draws_ab[1:n])
}
