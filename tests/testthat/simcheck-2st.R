
d <- .2; d_sd <- 1; n1 <- n2 <- 91; alpha <- .05
df <- n1 + n2 - 2
cv <- stats::qt(1 - alpha / 2, df)
d_draws <- rtruncnorm(1e6, a = 0, b = Inf, mean = d, sd = d_sd)
ncp_draws <- d_draws * sqrt(n1 * n2 / (n1 + n2))
pow_draws <- stats::pt(cv, df = df, ncp = ncp_draws, lower.tail = FALSE) +
  stats::pt(-cv, df = df, ncp = ncp_draws, lower.tail = TRUE)
# expected power
mean(pow_draws)
ep_2st(d_est = d, d_sd = d_sd,
       n1 = n1, n2 = n2, prior_d = "trunc_norm", trunc_d = c(0, Inf))
# assurance level of power
mean(pow_draws > .8)
al_2st(d_est = d, d_sd = d_sd,
       n1 = n1, n2 = n2, prior_d = "trunc_norm", trunc_d = c(0, Inf))

# solve for n
n_2st(d_est = .2, d_sd = 1, ep = .8,
      prior_d = "trunc_norm", trunc_d = c(0, Inf),
      test = "two.sided")
n_2st(d_est = .2, d_sd = 1, al = .8,
      prior_d = "trunc_norm", trunc_d = c(0, Inf),
      test = "two.sided")

#### additional examples

d <- .2; d_sd <- sqrt(1.144598 + 0.22837262); n1 <- n2 <- 300; alpha <- .05
df <- n1 + n2 - 2
cv <- stats::qt(1 - alpha / 2, df)
d_draws <- rtruncnorm(1e6, a = 0, b = Inf, mean = d, sd = d_sd)
ncp_draws <- d_draws * sqrt(n1 * n2 / (n1 + n2))
pow_draws <- stats::pt(cv, df = df, ncp = ncp_draws, lower.tail = FALSE) +
  stats::pt(-cv, df = df, ncp = ncp_draws, lower.tail = TRUE)
# expected power
mean(pow_draws)
ep_2st(d_est = d, d_sd = sqrt(1.144598 + 0.22837262),
       n1 = n1, n2 = n2, prior_d = "trunc_norm", trunc_d = c(0, Inf))
# assurance level of power
mean(pow_draws > .8)
al_2st(d_est = d, d_sd = sqrt(1.144598 + 0.22837262),
       n1 = n1, n2 = n2, prior_d = "trunc_norm", trunc_d = c(0, Inf))


# zero-inflated
ep_2st(d_est = .2, d_sd = sqrt(1.144598 + 0.22837262),
       n1 = 300, n2 = 300, prior_d = "zero-inflated", ndraws = 1e6)
n_2st(d_est = .2, d_sd = sqrt(0.22837262),
      prior_d = "zero-inflated", ndraws = 1e5)



n_2st(0.7579454, 0.2343928, al = .8,
      prior_d = "trunc_norm", trunc_d = c(0.4315931, Inf),
      test = "one.sided")
