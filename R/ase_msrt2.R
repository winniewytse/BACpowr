ase_msrt2 <- function(rho, rho_sd, omega, omega_sd,
                      J, n, se = 0.05, rsq1 = 0, rsq2 = 0,
                      K = 0, P = .5, ...) {

  if (rho_sd == 0) {
    if (omega_sd == 0) { # (1) rho_sd = omega_sd = 0
      se_msrt2(rho = rho, omega = omega, J = J, n = n,
               rsq1 = rsq1, rsq2 = rsq2, K = K, P = P)
    } else { # (2) rho_sd = 0
      omega_ab <- gamma_ab(omega, omega_sd)
      stats::pgamma(
        inv_se_msrt2(se = se, J = J, n = n, rho = rho,
                     rsq1 = rsq1, rsq2 = rsq2, K = K, P = P),
        shape = omega_ab[1], rate = omega_ab[2]
      )
    }
  } else {
    if (omega_sd == 0) { # (3) omega_sd = 0
      rho_ab <- get_ab(rho, rho_sd)
      stats::pbeta(
        inv_se_msrt2(se = se, J = J, n = n, omega = omega,
                     rsq1 = rsq1, rsq2 = rsq2, K = K, P = P),
        shape1 = rho_ab[1], shape2 = rho_ab[2]# , lower.tail = FALSE
      )
    } else { # (4)
      rho_ab <- get_ab(rho, rho_sd)
      omega_ab <- gamma_ab(omega, omega_sd)
      cubature::cuhre(
        function(x) {
          stats::pbeta(
            inv_se_msrt2(se = se, J = J, n = n, omega = x,
                         rsq1 = rsq1, rsq2 = rsq2, K = K, P = P),
            shape1 = rho_ab[1], shape2 = rho_ab[2]
          ) * stats::dgamma(x, shape = omega_ab[1], rate = omega_ab[2])
        },
        lowerLimit = 0, upperLimit = 1
      )$integral
    }
  }
}
