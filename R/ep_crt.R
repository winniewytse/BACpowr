ep_crt <- function(J, n, d_est, d_se = 0, 
                       rho_est, rho_se = 0, 
                       r2_est = 0, r2_se = 0, 
                       test = "two-tailed", 
                       relTol = 1e-3, absTol= 1e-50) {
  if (d_se < .005) {d_se = 0} else {d_se = d_se}
  if (d_se == 0) {
    if (rho_se == 0) {
      if (r2_se == 0) {             # (1) d_se = rho_se = r2_se = 0
        power_crt(J, n, rho_est, r2_est, d_est, test)
      } else {                      # (2) d_se = rho_se = 0
        r2_ab <- get_ab(r2_est, r2_se)
        cubature::cuhre(
          function(arg) {
            r2 <- arg[1]
            power_crt(J, n, rho_est, r2, d_est, test) * 
              dbeta(r2, r2_ab[1], r2_ab[2])
          }, 
          lowerLimit = 0, upperLimit = 1, 
          relTol = relTol, absTol= absTol
        )$integral
      }
    } else {                        
      if (r2_se == 0) {             # (3) d_se = r2_se = 0
        rho_ab <- get_ab(rho_est, rho_se)
        cubature::cuhre(
          function(arg) {
            rho <- arg[1]
            power_crt(J, n, rho, r2_est, d_est, test) * 
              dbeta(rho, rho_ab[1], rho_ab[2])
          }, 
          lowerLimit = 0, upperLimit = 1, 
          relTol = relTol, absTol= absTol
        )$integral
      } else {                      # (4) d_se = 0
        rho_ab <- get_ab(rho_est, rho_se)
        r2_ab <- get_ab(r2_est, r2_se)
        cubature::cuhre(
          function(arg) {
            rho <- arg[1]
            r2 <- arg[2]
            power_crt(J, n, rho, r2, d_est, test) *
              dbeta(rho_ab[1], rho_ab[2]) *
              dbeta(r2_ab[1], r2_ab[2])
          }, 
          lowerLimit = c(0, 0), upperLimit = c(1, 1), 
          relTol = relTol, absTol= absTol
        )$integral
      }
    }
  } else {                          
    if (rho_se == 0) {
      if (r2_se == 0) {             # (5) rho_se = r2_se = 0
        cubature::cuhre(
          function(arg) {
            delta <- arg[1]
            power_crt(J, n, rho_est, r2_est, delta, test) * 
              dnorm(delta, d_est, d_se)
          }, 
          lowerLimit = -Inf, upperLimit = Inf, 
          relTol = relTol, absTol= absTol
        )$integral
      } else {                      # (6) rho_se = 0
        r2_ab <- get_ab(r2_est, r2_se)
        cubature::cuhre(
          function(arg) {
            delta <- arg[1]
            r2 <- arg[2]
            power_crt(J, n, rho_est, r2, delta, test) * 
              dbeta(r2, r2_ab[1], r2_ab[2]) * 
              dnorm(delta, d_est, d_se)
          }, 
          lowerLimit = c(-Inf, 0), upperLimit = c(Inf, 1), 
          relTol = relTol, absTol= absTol
        )$integral
      }
    } else {
      if (r2_se == 0) {             # (7) r2_se = 0
        rho_ab <- get_ab(rho_est, rho_se)
        cubature::cuhre(
          function(arg) {
            delta <- arg[1]
            rho <- arg[2]
            power_crt(J, n, rho, r2_est, delta, test) * 
              dbeta(rho, rho_ab[1], rho_ab[2]) * 
              dnorm(delta, d_est, d_se)
          }, 
          lowerLimit = c(-Inf, 0), upperLimit = c(Inf, 1), 
          relTol = relTol, absTol= absTol
        )$integral
      } else {                      # (8) 
        rho_ab <- get_ab(rho_est, rho_se)
        r2_ab <- get_ab(r2_est, r2_se)
        cubature::cuhre(
          function(arg) {
            delta <- arg[1]
            rho <- arg[2]
            r2 <- arg[3]
            power_crt(J, n, rho, r2, delta, test) *
              dbeta(rho, rho_ab[1], rho_ab[2]) *
              dbeta(r2, r2_ab[1], r2_ab[2]) *
              dnorm(delta, d_est, d_se)
          }, 
          lowerLimit = c(-Inf, 0, 0), upperLimit = c(Inf, 1, 1), 
          relTol = relTol, absTol= absTol
        )$integral
      }
    }
  }
}