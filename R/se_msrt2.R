se_msrt2 <- function(rho, omega, J, n, rsq1 = 0, rsq2 = 0,
                     K = 0, P = .5) {
  sqrt((rho * omega * (1 - rsq2) * P * (1 - P) * n +
          (1 - rho) * (1 - rsq1)) /
         (P * (1 - P) * J * n))
}
