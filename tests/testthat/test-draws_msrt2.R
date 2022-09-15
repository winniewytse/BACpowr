
draws_d <- rnorm(8000, mean = .2, sd = .1)
rho_ab <- get_ab(.1, .05)
draws_rho <- rbeta(8000, shape1 = rho_ab[1], shape2 = rho_ab[2])
omega_ab <- gamma_ab(.3, .1)
draws_omega <- rgamma(8000, shape = omega_ab[1], omega_ab[2])

draws_msrt2(draws_d = draws_d, draws_rho = draws_rho, draws_omega = draws_omega,
            J = 75, n = 20, power = .8, goal = "ep")
draws_msrt2(draws_d = draws_d, draws_rho = draws_rho, draws_omega = draws_omega,
            J = 126, n = 20, power = .8, goal = "al")
draws_msrt2(draws_d = draws_d, draws_rho = draws_rho, draws_omega = draws_omega,
            J = 41, n = 20, precision = .15, goal = "apr")


Jn_msrt2_draws(draws_d = draws_d, draws_rho = draws_rho, draws_omega = draws_omega,
               n = 20, power = .8, ep = .8)
Jn_msrt2_draws(draws_d = draws_d, draws_rho = draws_rho, draws_omega = draws_omega,
               n = 20, power = .8, al = .8)
Jn_msrt2_draws(draws_d = draws_d, draws_rho = draws_rho, draws_omega = draws_omega,
               n = 20, precision = .15, apr = .8)
