
# uncertainty in rho
al_crt2(120, 20, .3, 0, .2, .1)
ep_crt2(120, 20, .3, 0, .2, .1)

df <- 120 - 2
shapes <- get_ab(.2, .1)
rho_draws <- rbeta(1e6, shapes[1], shapes[2])
ncp_draws <- .3 * sqrt(120 * 20 / 4 / (1 + (20 - 1) * rho_draws))
pow_draws <- pt(qt(.975, df), df = df, ncp = ncp_draws, lower.tail = FALSE) +
  pt(-qt(.975, df), df = df, ncp = ncp_draws, lower.tail = TRUE)
mean(pow_draws > .8)
mean(pow_draws)


# uncertainty in delta
al_crt2(120, 20, d_est = .3, d_sd = .1, rho_est = .2, rho_sd = 0)
ep_crt2(120, 20, .3, .1, .2, 0)

df <- 120 - 2
shapes <- get_ab(.2, .1)
delta_draws <- rnorm(1e6, .3, .1)
ncp_draws <- delta_draws * sqrt(120 * 20 / 4 / (1 + (20 - 1) * .2))
pow_draws <- pt(qt(.975, df), df = df, ncp = ncp_draws, lower.tail = FALSE) +
  pt(-qt(.975, df), df = df, ncp = ncp_draws, lower.tail = TRUE)
mean(pow_draws > .8)
mean(pow_draws)


# uncertainty in both rho and delta
al_crt2(120, 20, .3, .8, .2, .1)
ep_crt2(120, 20, .3, .8, .2, .1)

df <- 120 - 2
shapes <- get_ab(.2, .1)
rho_draws <- rbeta(1e6, shapes[1], shapes[2])
delta_draws <- rnorm(1e6, .3, .8)
ncp_draws <- delta_draws * sqrt(120 * 20 / 4 / (1 + (20 - 1) * rho_draws))
pow_draws <- pt(qt(.975, df), df = df, ncp = ncp_draws, lower.tail = FALSE) +
  pt(-qt(.975, df), df = df, ncp = ncp_draws, lower.tail = TRUE)
mean(pow_draws > .8)
mean(pow_draws)


# additional check
al_crt2(120, 20, .3, .1, .2, .2)
ep_crt2(120, 20, .3, .1, .2, .2)

df <- 120 - 2
shapes <- get_ab(.2, .2)
rho_draws <- rbeta(5e6, shapes[1], shapes[2])
delta_draws <- rnorm(5e6, .3, .1)
ncp_draws <- delta_draws * sqrt(120 * 20 / 4 / (1 + (20 - 1) * rho_draws))
pow_draws <- pt(qt(.975, df), df = df, ncp = ncp_draws, lower.tail = FALSE) +
  pt(-qt(.975, df), df = df, ncp = ncp_draws, lower.tail = TRUE)
mean(pow_draws > .8)
mean(pow_draws)
