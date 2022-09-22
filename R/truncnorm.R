#' Truncated normal distribution
#' @export

dtruncnorm <- function(x, a = -Inf, b = Inf, mean = 0, sd = 1, ...) {
  if (x <= a | x >= b) {
    0
  } else {
    dnorm(x, mean = mean, sd = sd, ...) /
      (pnorm(b, mean = mean, sd = sd, ...) - pnorm(a, mean = mean, sd = sd, ...))
  }
}

ptruncnorm <- function(q, a = -Inf, b = Inf, mean = 0, sd = 1, ...) {
  if (q <= a) 0
  else if (q >= b) 1
  else
    (pnorm(q, mean, sd, ...) - pnorm(a, mean, sd, ...)) /
    (pnorm(b, mean, sd, ...) - pnorm(a, mean, sd, ...))
}

qtruncnorm <- function(p, a = -Inf, b = Inf, mean = 0, sd = 1, ...) {
  qnorm(pnorm(a, mean, sd, ...) + p *
          (pnorm(b, mean, sd, ...) - pnorm(a, mean, sd, ...)),
        mean, sd, ...)
}

rtruncnorm <- function(n, a = -Inf, b = Inf, mean = 0, sd = 1, ...) {
  ndraws <- 0
  draws_ab <- c()
  while (ndraws < n) {
    draws <- rnorm(n, mean, sd, ...)
    draws_ab <- c(draws_ab, draws[draws > a & draws < b])
    ndraws <- length(draws_ab)
  }
  return(draws_ab[1:n])
}
