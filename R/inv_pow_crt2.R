inv_pow_crt2 <- function(power, J, n, d_est = NULL, rho_est = NULL,
                         r2_est = NULL, K = 0) {
  df <- J - K - 2
  cv <- stats::qt(.975, df)
  # cv <- qnorm(.975)
  if (is.null(d_est)) {
    (cv - qnorm(1 - power)) *
      sqrt((4 * (1 + (n * (1 - r2_est) - 1) * rho_est)) / (J * n))
  } else if (is.null(rho_est)) {
    (J * n / 4 * d_est^2 / (cv - qnorm(1 - power))^2 - 1) / (n * (1 - r2_est) - 1)
  } else if (is.null(r2_est)) {
    1 - ((J * n / 4 * d_est^2 / (cv - qnorm(1 - power))^2 - 1) / rho_est + 1) / n
  }
}
