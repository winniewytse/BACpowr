ese_msrt2 <- function(rho, rho_sd, omega, omega_sd,
                      J, n, se = 0.05, rsq1 = 0, rsq2 = 0,
                      K = 0, P = .5, ...) {
  df <- J - K - 1
  if (rho_sd == 0) {
    if (omega_sd == 0) {
      se_msrt2(rho = rho, omega = omega, J = J, n = n,
               rsq1 = rsq1, rsq2 = rsq2, K = K, P = P)
    } else {
      omega_ab <- gamma_ab(omega, omega_sd)
      cubature::hcubature(
        function(omega_i) {
          se_msrt2(rho = rho, omega = omega_i, J = J, n = n,
                   rsq1 = rsq1, rsq2 = rsq2, K = K, P = P) *
            stats::dgamma(omega_i, omega_ab[1], omega_ab[2])
        },
        lowerLimit = 0, upperLimit = 1,
        vectorInterface = TRUE
      )$integral
    }
  } else {
    if (omega_sd == 0) {
      rho_ab <- get_ab(rho, rho_sd)
      cubature::hcubature(
        function(rho_i) {
          se_msrt2(rho = rho_i, omega = omega, J = J, n = n,
                   rsq1 = rsq1, rsq2 = rsq2, K = K, P = P) *
            stats::dbeta(rho_i, rho_ab[1], rho_ab[2])
        },
        lowerLimit = 0, upperLimit = 1,
        vectorInterface = TRUE, ...
      )$integral
    } else {
      rho_ab <- get_ab(rho, rho_sd)
      omega_ab <- gamma_ab(omega, omega_sd)
      cubature::hcubature(
        function(x) {
          matrix(apply(x, 2, function(args) {
            rho_i <- args[1]
            omega_i <- args[2]
            se_msrt2(rho = rho_i, omega = omega_i, J = J, n = n,
                     rsq1 = rsq1, rsq2 = rsq2, K = K, P = P) *
              stats::dbeta(rho_i, rho_ab[1], rho_ab[2]) *
              stats::dgamma(omega_i, omega_ab[1], omega_ab[2])
          }), ncol = ncol(x))
        },
        lowerLimit = c(0, 0), upperLimit = c(1, 1),
        vectorInterface = TRUE, ...
      )$integral
    }
  }
}
