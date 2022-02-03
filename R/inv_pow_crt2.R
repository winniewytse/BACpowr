inv_pow_crt2 <- function(power, J, n, d_est = NULL, rho_est = NULL,
                         r2_est = NULL, K = 0, test = "two-tailed") {
  if (test == "two-tailed") {
    if (is.null(J)) stop("Missing J")
    if (is.null(n)) stop("Missing n")
    df1 <- 1
    df2 <- J - K - 2
    cv <- stats::qf(.95, df1, df2)
    if (is.null(d_est)) {
      inv <- function(d_est) {
        ncp <- d_est^2 * (J * n / 4 / (1 + (n * (1 - r2_est) - 1) * rho_est))
        pf(cv, df1, df2, ncp, lower.tail = FALSE) - power
      }
      stats::uniroot(inv, c(0, 10))$root
    } else if (is.null(rho_est)) {
      inv <- function(rho_est) {
        ncp <- d_est^2 * (J * n / 4 / (1 + (n * (1 - r2_est) - 1) * rho_est))
        pf(cv, df1, df2, ncp, lower.tail = FALSE) - power
      }
      stats::uniroot(inv, c(0, 1))$root
    } else if (is.null(r2_est)) {
      inv <- function(r2_est) {
        ncp <- d_est^2 * (J * n / 4 / (1 + (n * (1 - r2_est) - 1) * rho_est))
        pf(cv, df1, df2, ncp, lower.tail = FALSE) - power
      }
      stats::uniroot(inv, c(0, 1))$root
    }
  } else if (test == "one-tailed") {
    df <- J - K - 2
    cv <- stats::qt(.95, df)
    if (is.null(d_est)) {
      (cv - qnorm(1 - power)) *
        sqrt((4 * (1 + (n * (1 - r2_est) - 1) * rho_est)) / (J * n))
    } else if (is.null(rho_est)) {
      (J * n / 4 * d_est^2 / (cv - qnorm(1 - power))^2 - 1) / (n * (1 - r2_est) - 1)
    } else if (is.null(r2_est)) {
      1 - ((J * n / 4 * d_est^2 / (cv - qnorm(1 - power))^2 - 1) / rho_est + 1) / n
    }
  }
}
