
#### Two-level CRT ####

# uncertainty in rho
al_crt2(J = 120, n = 20, d_est = .3, d_sd = 0, rho_est = .2, rho_sd = .1)
ep_crt2(J = 120, n = 20, d_est = .3, d_sd = 0, rho_est = .2, rho_sd = .1)

df <- 120 - 2
shapes <- get_ab(.2, .1)
rho_draws <- rbeta(1e6, shapes[1], shapes[2])
ncp_draws <- .3 * sqrt(120 * 20 / 4 / (1 + (20 - 1) * rho_draws))
pow_draws <- pt(qt(.975, df), df = df, ncp = ncp_draws, lower.tail = FALSE) +
  pt(-qt(.975, df), df = df, ncp = ncp_draws, lower.tail = TRUE)
mean(pow_draws > .8)
mean(pow_draws)


# uncertainty in delta (two-sided)
al_crt2(J = 120, n = 20, d_est = .3, d_sd = .1, rho_est = .2, rho_sd = 0)
ep_crt2(J = 120, n = 20, d_est = .3, d_sd = .1, rho_est = .2, rho_sd = 0)

df <- 120 - 2
shapes <- get_ab(.2, .1)
delta_draws <- rnorm(1e6, .3, .1)
ncp_draws <- delta_draws * sqrt(120 * 20 / 4 / (1 + (20 - 1) * .2))
pow_draws <- pt(qt(.975, df), df = df, ncp = ncp_draws, lower.tail = FALSE) +
  pt(-qt(.975, df), df = df, ncp = ncp_draws, lower.tail = TRUE)
mean(pow_draws > .8)
mean(pow_draws)

# uncertainty in delta (one-sided)
al_crt2(J = 120, n = 20, d_est = .3, d_sd = .1, rho_est = .2, rho_sd = 0, test = "one.sided")
ep_crt2(J = 120, n = 20, d_est = .3, d_sd = .1, rho_est = .2, rho_sd = 0)

df <- 120 - 2
shapes <- get_ab(.2, .1)
delta_draws <- rnorm(1e6, .3, .1)
ncp_draws <- delta_draws * sqrt(120 * 20 / 4 / (1 + (20 - 1) * .2))
pow_draws <- pt(qt(.95, df), df = df, ncp = ncp_draws, lower.tail = FALSE)
mean(pow_draws > .8)
mean(pow_draws)


# uncertainty in both rho and delta
al_crt2(J = 120, n = 20, d_est = .3, d_sd = .8, rho_est = .2, rho_sd = .1)
ep_crt2(J = 120, n = 20, d_est = .3, d_sd = .8, rho_est = .2, rho_sd = .1)

J <- 120; n <- 20
d <- .3; d_sd <- .8; rho <- .2; rho_sd <- .1
df <- J - 2
shapes <- get_ab(rho, rho_sd)
rho_draws <- rbeta(1e6, shapes[1], shapes[2])
delta_draws <- rnorm(1e6, d, d_sd)
ncp_draws <- delta_draws * sqrt(J * n / 4 / (1 + (n - 1) * rho_draws))
pow_draws <- pt(qt(.975, df), df = df, ncp = ncp_draws, lower.tail = FALSE) +
  pt(-qt(.975, df), df = df, ncp = ncp_draws, lower.tail = TRUE)
mean(pow_draws > .8)
mean(pow_draws)


# additional check
al_crt2(J = 120, n = 20, d_est = .3, d_sd = .1, rho_est = .2, rho_sd = .2)
ep_crt2(J = 120, n = 20, d_est = .3, d_sd = .1, rho_est = .2, rho_sd = .2)

df <- 120 - 2
shapes <- get_ab(.2, .2)
rho_draws <- rbeta(5e6, shapes[1], shapes[2])
delta_draws <- rnorm(5e6, .3, .1)
ncp_draws <- delta_draws * sqrt(120 * 20 / 4 / (1 + (20 - 1) * rho_draws))
pow_draws <- pt(qt(.975, df), df = df, ncp = ncp_draws, lower.tail = FALSE) +
  pt(-qt(.975, df), df = df, ncp = ncp_draws, lower.tail = TRUE)
mean(pow_draws > .8)
mean(pow_draws)


#### Two-level MSRT ####

# uncertainty in delta
al_msrt2(J = 50, n = 30, d_est = .3, d_sd = .1, rho_est = .2, rho_sd = 0,
         omega_est = .3, omega_sd = 0)
ep_msrt2(J = 50, n = 30, d_est = .3, d_sd = .1, rho_est = .2, rho_sd = 0,
         omega_est = .3, omega_sd = 0)

J <- 50; n <- 30; P <- .5; rsq1 <- 0; rsq2 <- 0
d <- .3; d_sd = .1; rho = .2; omega <- .3
df <- J - 1
d_draws <- rnorm(1e6, d, d_sd)
ncp_draws <- d_draws * sqrt(P * (1 - P) * J * n /
                              (rho * omega * (1 - rsq2) * P * (1 - P) * n +
                                 (1 - rho) * (1 - rsq1)))
pow_draws <- pt(qt(.975, df), df = df, ncp = ncp_draws, lower.tail = FALSE) +
  pt(-qt(.975, df), df = df, ncp = ncp_draws, lower.tail = TRUE)
mean(pow_draws > .8)
mean(pow_draws)

# uncertainty in rho
al_msrt2(J = 20, n = 25, d_est = .3, d_sd = 0, rho_est = .2, rho_sd = .1,
         omega_est = .3, omega_sd = 0)
ep_msrt2(J = 20, n = 25, d_est = .3, d_sd = 0, rho_est = .2, rho_sd = .1,
         omega_est = .3, omega_sd = 0)

J <- 20; n <- 25; P <- .5; rsq1 <- 0; rsq2 <- 0
d <- .3; rho <- .2; rho_sd <- .1; omega <- .3
df <- J - 1
rho_ab <- get_ab(rho, rho_sd)
rho_draws <- rbeta(1e6, rho_ab[1], rho_ab[2])
ncp_draws <- d * sqrt(P * (1 - P) * J * n /
                        (rho_draws * omega * (1 - rsq2) * P * (1 - P) * n +
                           (1 - rho_draws) * (1 - rsq1)))
pow_draws <- pt(qt(.975, df), df = df, ncp = ncp_draws, lower.tail = FALSE) +
  pt(-qt(.975, df), df = df, ncp = ncp_draws, lower.tail = TRUE)
mean(pow_draws > .8)
mean(pow_draws)


# uncertainty in omega
al_msrt2(J = 20, n = 23, d_est = .3, d_sd = 0, rho_est = .2, rho_sd = 0,
         omega_est = .3, omega_sd = .1)
ep_msrt2(J = 20, n = 23, d_est = .3, d_sd = 0, rho_est = .2, rho_sd = 0,
         omega_est = .3, omega_sd = .1)

J <- 20; n <- 23; P <- .5; rsq1 <- 0; rsq2 <- 0
d <- .3; rho <- .2; omega <- .3; omega_sd <- .1
df <- J - 1
omega_ab <- gamma_ab(omega, omega_sd)
omega_draws <- rgamma(1e6, omega_ab[1], omega_ab[2])
ncp_draws <- d * sqrt(P * (1 - P) * J * n /
                        (rho * omega_draws * (1 - rsq2) * P * (1 - P) * n +
                           (1 - rho) * (1 - rsq1)))
pow_draws <- pt(qt(.975, df), df = df, ncp = ncp_draws, lower.tail = FALSE) +
  pt(-qt(.975, df), df = df, ncp = ncp_draws, lower.tail = TRUE)
mean(pow_draws > .8)
mean(pow_draws)

# uncertainty in delta and rho
al_msrt2(J = 50, n = 30, d_est = .3, d_sd = .1, rho_est = .2, rho_sd = .1,
         omega_est = .3, omega_sd = 0)
ep_msrt2(J = 50, n = 30, d_est = .3, d_sd = .1, rho_est = .2, rho_sd = .1,
         omega_est = .3, omega_sd = 0)

J <- 50; n <- 30; P <- .5; rsq1 <- 0; rsq2 <- 0
d <- .3; d_sd = .1; rho = .2; rho_sd = .1; omega <- .3; omega_sd <- 0
df <- J - 1
d_draws <- rnorm(1e6, d, d_sd)
rho_ab <- get_ab(rho, rho_sd)
rho_draws <- rbeta(1e6, rho_ab[1], rho_ab[2])
ncp_draws <- d_draws * sqrt(P * (1 - P) * J * n /
                              (rho_draws * omega * (1 - rsq2) * P * (1 - P) * n +
                                 (1 - rho_draws) * (1 - rsq1)))
pow_draws <- pt(qt(.975, df), df = df, ncp = ncp_draws, lower.tail = FALSE) +
  pt(-qt(.975, df), df = df, ncp = ncp_draws, lower.tail = TRUE)
mean(pow_draws > .8)
mean(pow_draws)

# uncertainty in delta and omega
al_msrt2(J = 50, n = 30, d_est = .3, d_sd = .1, rho_est = .2, rho_sd = 0,
         omega_est = .3, omega_sd = .1)
ep_msrt2(J = 50, n = 30, d_est = .3, d_sd = .1, rho_est = .2, rho_sd = 0,
         omega_est = .3, omega_sd = .1)

J <- 50; n <- 30; P <- .5; rsq1 <- 0; rsq2 <- 0
d <- .3; d_sd = .1; rho = .2; omega <- .3; omega_sd <- .1
df <- J - 1
d_draws <- rnorm(1e6, d, d_sd)
omega_ab <- gamma_ab(omega, omega_sd)
omega_draws <- rgamma(1e6, omega_ab[1], omega_ab[2])
ncp_draws <- d_draws * sqrt(P * (1 - P) * J * n /
                              (rho * omega_draws * (1 - rsq2) * P * (1 - P) * n +
                                 (1 - rho) * (1 - rsq1)))
pow_draws <- pt(qt(.975, df), df = df, ncp = ncp_draws, lower.tail = FALSE) +
  pt(-qt(.975, df), df = df, ncp = ncp_draws, lower.tail = TRUE)
mean(pow_draws > .8)
mean(pow_draws)

# uncertainty in rho and omega
al_msrt2(J = 20, n = 50, d_est = .3, d_sd = 0, rho_est = .3, rho_sd = .1,
         omega_est = .5, omega_sd = .1)
ep_msrt2(J = 20, n = 50, d_est = .3, d_sd = 0, rho_est = .3, rho_sd = .1,
         omega_est = .5, omega_sd = .1)

J <- 20; n <- 50; P <- .5; rsq1 <- 0; rsq2 <- 0
d <- .3; rho = .3; rho_sd <- .1; omega <- .5; omega_sd <- .1
df <- J - 1
rho_ab <- get_ab(rho, rho_sd)
rho_draws <- rbeta(1e6, rho_ab[1], rho_ab[2])
omega_ab <- gamma_ab(omega, omega_sd)
omega_draws <- rgamma(1e6, omega_ab[1], omega_ab[2])
ncp_draws <- d * sqrt(P * (1 - P) * J * n /
                              (rho_draws * omega_draws * (1 - rsq2) * P * (1 - P) * n +
                                 (1 - rho_draws) * (1 - rsq1)))
pow_draws <- pt(qt(.975, df), df = df, ncp = ncp_draws, lower.tail = FALSE) +
  pt(-qt(.975, df), df = df, ncp = ncp_draws, lower.tail = TRUE)
mean(pow_draws > .8)
mean(pow_draws)


# uncertainty in delta, rho, and omega
al_msrt2(J = 20, n = 50, d_est = .3, d_sd = .1, rho_est = .3, rho_sd = .1,
         omega_est = .5, omega_sd = .1)
ep_msrt2(J = 20, n = 50, d_est = .3, d_sd = .1, rho_est = .3, rho_sd = .1,
         omega_est = .5, omega_sd = .1)

J <- 20; n <- 50; P <- .5; rsq1 <- 0; rsq2 <- 0
d <- .3; d_sd <- .1; rho = .3; rho_sd <- .1; omega <- .5; omega_sd <- .1
df <- J - 1
d_draws <- rnorm(1e6, d, d_sd)
rho_ab <- get_ab(rho, rho_sd)
rho_draws <- rbeta(1e6, rho_ab[1], rho_ab[2])
omega_ab <- gamma_ab(omega, omega_sd)
omega_draws <- rgamma(1e6, omega_ab[1], omega_ab[2])
ncp_draws <- d_draws * sqrt(P * (1 - P) * J * n /
                        (rho_draws * omega_draws * (1 - rsq2) * P * (1 - P) * n +
                           (1 - rho_draws) * (1 - rsq1)))
pow_draws <- pt(qt(.975, df), df = df, ncp = ncp_draws, lower.tail = FALSE) +
  pt(-qt(.975, df), df = df, ncp = ncp_draws, lower.tail = TRUE)
mean(pow_draws > .8)
mean(pow_draws)


#### Independent Sample t-test ####

al_2st(231, 231, .6360188, sqrt(0.19873))

n1 <- n2 <- 231
d <- .6360188; d_sd <- sqrt(0.19873)
df <- n1 + n2 - 2
d_draws <- rnorm(1e6, d, d_sd)
ncp_draws <- d_draws * sqrt(n1 * n2 / (n1 + n2))
pow_draws <- pt(qt(.975, df), df = df, ncp = ncp_draws, lower.tail = FALSE) +
  pt(-qt(.975, df), df = df, ncp = ncp_draws, lower.tail = TRUE)
mean(pow_draws > .8)
mean(pow_draws)
