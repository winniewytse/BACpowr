
# uncertainty in rho
ese_msrt2(rho = .2, rho_sd = .1, omega = .3, omega_sd = 0, J = 50, n = 30)
ase_msrt2(rho = .2, rho_sd = .1, omega = .3, omega_sd = 0, se = .06, J = 50, n = 30)

J <- 50; n <- 30; P <- .5; rsq1 <- 0; rsq2 <- 0; K <- 0
rho <- .2; rho_sd <- .1; omega <- .3
rho_ab <- get_ab(rho, rho_sd)
rho_draws <- rbeta(1e6, rho_ab[1], rho_ab[2])
se_draws <- sqrt((rho_draws * omega * (1 - rsq2) * P * (1 - P) * n +
                    (1 - rho_draws) * (1 - rsq1)) / (P * (1 - P) * J * n))
mean(se_draws)
mean(se_draws < .06)


# uncertainty in omega
ese_msrt2(rho = .2, rho_sd = 0, omega = .3, omega_sd = .1, J = 50, n = 30)
ase_msrt2(rho = .2, rho_sd = 0, omega = .3, omega_sd = .1, se = .06, J = 50, n = 30)

J <- 50; n <- 30; P <- .5; rsq1 <- 0; rsq2 <- 0; K <- 0
rho <- .2; omega <- .3; omega_sd <- .1
omega_ab <- gamma_ab(omega, omega_sd)
omega_draws <- rgamma(1e6, omega_ab[1], omega_ab[2])
se_draws <- sqrt((rho * omega_draws * (1 - rsq2) * P * (1 - P) * n +
                    (1 - rho) * (1 - rsq1)) / (P * (1 - P) * J * n))
mean(se_draws)
mean(se_draws < .06)



# uncertainty in rho and omega
ese_msrt2(rho = .2, rho_sd = .1, omega = .3, omega_sd = .1, J = 50, n = 30)
ase_msrt2(rho = .2, rho_sd = .1, omega = .3, omega_sd = .1, se = .06, J = 50, n = 30)

J <- 50; n <- 30; P <- .5; rsq1 <- 0; rsq2 <- 0; K <- 0
rho <- .2; rho_sd <- .1; omega <- .3; omega_sd <- .1
rho_ab <- get_ab(rho, rho_sd)
rho_draws <- rbeta(1e6, rho_ab[1], rho_ab[2])
omega_ab <- gamma_ab(omega, omega_sd)
omega_draws <- rgamma(1e6, omega_ab[1], omega_ab[2])
se_draws <- sqrt((rho_draws * omega_draws * (1 - rsq2) * P * (1 - P) * n +
                    (1 - rho_draws) * (1 - rsq1)) / (P * (1 - P) * J * n))
mean(se_draws)
mean(se_draws < .06)


