
# uncertainty in rho
# ese_msrt2(rho = .2, rho_sd = .1, omega = .3, omega_sd = 0, J = 50, n = 30)
apr_msrt2(rho = .2, rho_sd = .1, omega = .3, omega_sd = 0, precision = .12,
          J = 50, n = 30)

J <- 50; n <- 30; P <- .5; rsq1 <- 0; rsq2 <- 0; K <- 0; alpha <- .05
rho <- .2; rho_sd <- .1; omega <- .3
rho_ab <- get_ab(rho, rho_sd)
rho_draws <- rbeta(1e6, rho_ab[1], rho_ab[2])
df <- J - K - 1
t_crit <- qt(1 - alpha / 2, df = df)
prec_draws <- t_crit * sqrt((rho_draws * omega * (1 - rsq2) * P * (1 - P) * n +
                               (1 - rho_draws) * (1 - rsq1)) / (P * (1 - P) * J * n))
# mean(prec_draws)
mean(prec_draws < .12)

apr_msrt2(J = 49.61714, n = 40, P = .4,
          rho = 0.4563287, rho_sd = .21,
          omega = 0.2153255, omega_sd = 0, precision = .15)

J <- 49.61714; n <- 40; P <- .4; rsq1 <- 0; rsq2 <- 0; K <- 0; alpha <- .05
rho <- 0.4563287; rho_sd <- .21; omega <- 0.2153255
rho_ab <- get_ab(rho, rho_sd)
rho_draws <- rbeta(1e6, rho_ab[1], rho_ab[2])
df <- J - K - 1
t_crit <- qt(1 - alpha / 2, df = df)
prec_draws <- t_crit * sqrt((rho_draws * omega * (1 - rsq2) * P * (1 - P) * n +
                               (1 - rho_draws) * (1 - rsq1)) / (P * (1 - P) * J * n))
# mean(prec_draws)
mean(prec_draws < .15)


# uncertainty in omega
# ese_msrt2(rho = .2, rho_sd = 0, omega = .3, omega_sd = .1, J = 50, n = 30)
apr_msrt2(rho = .2, rho_sd = 0, omega = .3, omega_sd = .1, precision = .12,
          J = 50, n = 30)

J <- 50; n <- 30; P <- .5; rsq1 <- 0; rsq2 <- 0; K <- 0
rho <- .2; omega <- .3; omega_sd <- .1
omega_ab <- gamma_ab(omega, omega_sd)
omega_draws <- rgamma(1e6, omega_ab[1], omega_ab[2])
df <- J - K - 1
t_crit <- qt(1 - alpha / 2, df = df)
prec_draws <- t_crit * sqrt((rho * omega_draws * (1 - rsq2) * P * (1 - P) * n +
                    (1 - rho) * (1 - rsq1)) / (P * (1 - P) * J * n))
# mean(prec_draws)
mean(prec_draws < .12)



# uncertainty in rho and omega
# ese_msrt2(rho = .2, rho_sd = .1, omega = .3, omega_sd = .1, J = 50, n = 30)
apr_msrt2(rho = .2, rho_sd = .1, omega = .3, omega_sd = .1, precision = .12,
          J = 50, n = 30)

J <- 50; n <- 30; P <- .5; rsq1 <- 0; rsq2 <- 0; K <- 0; alpha <- .05
rho <- .2; rho_sd <- .1; omega <- .3; omega_sd <- .1
rho_ab <- get_ab(rho, rho_sd)
rho_draws <- rbeta(1e6, rho_ab[1], rho_ab[2])
omega_ab <- gamma_ab(omega, omega_sd)
omega_draws <- rgamma(1e6, omega_ab[1], omega_ab[2])
omega_draws2 <- omega_draws[omega_draws <= 1]
df <- J - K - 1
t_crit <- qt(1 - alpha / 2, df = df)
prec_draws <- t_crit * sqrt((rho_draws * omega_draws * (1 - rsq2) * P * (1 - P) * n +
                    (1 - rho_draws) * (1 - rsq1)) / (P * (1 - P) * J * n))
# mean(se_draws)
mean(prec_draws < .12)


apr_msrt2(rho = 0.4563287, rho_sd = 0.05040701, omega = 0.2153255,
          omega_sd = 0.03668986, precision = 0.15, J = 41, n = 20)


